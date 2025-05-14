#! /usr/bin/env python3

"""
Check a GFF3 file for errors and unwanted features, with an option to correct the errors and output a valid GFF3 file.

Count the number of Ns in each feature, remove features with N count greater than the specified threshold. (Requires FASTA)
Check and remove features with an end coordinates larger than the landmark sequence length. (Requires FASTA or ##sequence-region)
Check if the ##sequence-region matches the FASTA file. (Requires FASTA and ##sequence-region)
Add the ##sequence-region directives if missing. (Requires FASTA)
Check and correct the phase for CDS features.
"""
from collections import defaultdict
from itertools import groupby
try:
    from urllib import quote, unquote
except ImportError:
    from urllib.parse import quote, unquote
import re
import gzip
import string
import logging

logger = logging.getLogger(__name__)
# log.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)



ERROR_INFO = {
    'Ema0001': 'Parent feature start and end coordinates exceed those of child features',
    'Ema0002': 'Protein sequence contains internal stop codons',
    'Ema0003': 'This feature is not contained within the parent feature coordinates',
    'Ema0004': 'Incomplete gene feature that should contain at least one mRNA, exon, and CDS',
    'Ema0005': 'Pseudogene has invalid child feature type',
    'Ema0006': 'Wrong phase',
    'Ema0007': 'CDS and parent feature on different strands',
    'Ema0008': 'Warning for distinct isoforms that do not share any regions',
    'Ema0009': 'Incorrectly merged gene parent? Isoforms that do not share coding sequences are found',
    'Emr0001': 'Duplicate transcript found', # Error message has to be modified in lib/inter_model/inter_model.py
    'Emr0002': 'Incorrectly split gene parent?',
    'Emr0003': 'Duplicate ID',
    'Esf0001': 'Feature type may need to be changed to pseudogene',
    'Esf0002': 'Start/Stop is not a valid 1-based integer coordinate',
    'Esf0003': 'Strand information missing',
    'Esf0004': 'Seqid not found in any ##sequence-region',
    'Esf0005': 'Start is less than the ##sequence-region start',
    'Esf0006': 'End is greater than the ##sequence-region end',
    'Esf0007': 'Seqid not found in the embedded ##FASTA',
    'Esf0008': 'End is greater than the embedded ##FASTA sequence length',
    #'Esf0009': 'Found Ns in a feature using the embedded ##FASTA', # This error would be detected by the program as well, but the error message has to be modified in gff3.py
    'Esf0010': 'Seqid not found in the external FASTA file',
    'Esf0011': 'End is greater than the external FASTA sequence length',
    'Esf0012': 'Found Ns in a feature using the external FASTA', # This error would be detected by the program as well, but the error message has to be modified in gff3.py
    'Esf0013': 'White chars not allowed at the start of a line',
    'Esf0014': '##gff-version" missing from the first line',
    'Esf0015': 'Expecting certain fields in the feature',
    'Esf0016': '##sequence-region seqid may only appear once',
    'Esf0017': 'Start/End is not a valid integer',
    'Esf0018': 'Start is not less than or equal to end',
    'Esf0019': 'Version is not "3"',
    'Esf0020': 'Version is not a valid integer',
    'Esf0021': 'Unknown directive',
    'Esf0022': 'Features should contain 9 fields',
    'Esf0023': 'Escape certain characters',
    'Esf0024': 'Score is not a valid floating point number',
    'Esf0025': 'Strand has illegal characters',
    'Esf0026': 'Phase is not 0, 1, or 2, or not a valid integer',
    'Esf0027': 'Phase is required for all CDS features',
    'Esf0028': 'Attributes must escape the percent (%) sign and any control characters',
    'Esf0029': 'Attributes must contain one and only one equal (=) sign',
    'Esf0030': 'Empty attribute tag',
    'Esf0031': 'Empty attribute value',
    'Esf0032': 'Found multiple attribute tags',
    'Esf0033': 'Found ", " in a attribute, possible unescaped',
    'Esf0034': 'Attribute has identical values (count, value)',
    'Esf0035': 'Attribute has unresolved forward reference',
    'Esf0036': 'Value of a attribute contains unescaped ","',
    'Esf0037': 'Target attribute should have 3 or 4 values',
    'Esf0038': 'Start/End value of Target attribute is not a valid integer coordinate',
    'Esf0039': 'Strand value of Target attribute has illegal characters',
    'Esf0040': 'Value of Is_circular attribute is not "true"',
    'Esf0041': 'Unknown reserved (uppercase) attribute'
}

IDrequired = ['gene', 'pseudogene', 'mRNA', 'pseudogenic_transcript']

try:
    COMPLEMENT_TRANS = string.maketrans('TAGCtagc', 'ATCGATCG')
except AttributeError:
    COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
def complement(seq):
    return seq.translate(COMPLEMENT_TRANS)

BASES = ['t', 'c', 'a', 'g']
CODONS = [a+b+c for a in BASES for b in BASES for c in BASES]
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))
def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''
    for i in range(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = CODON_TABLE.get(codon, '!')
        if amino_acid != '!': # end of seq
            peptide += amino_acid
    return peptide

def fasta_file_to_dict(fasta_file, id=True, header=False, seq=False):
    """Returns a dict from a fasta file and the number of sequences as the second return value.
    fasta_file can be a string path or a file object.
    The key of fasta_dict can be set using the keyword arguments and
    results in a combination of id, header, sequence, in that order. joined with '||'. (default: id)
    Duplicate keys are checked and a warning is logged if found.
    The value of fasta_dict is a python dict with 3 keys: header, id and seq

    Changelog:
    2014/11/17:
    * Added support for url escaped id
    """
    fasta_file_f = fasta_file
    if isinstance(fasta_file, str):
        if fasta_file.endswith('.gz'):
            fasta_file_f = gzip.open(fasta_file, 'rt')  # gzip support
        else:
            fasta_file_f = open(fasta_file, 'r')

    fasta_dict = OrderedDict()
    keys = ['id', 'header', 'seq']
    flags = dict([('id', id), ('header', header), ('seq', seq)])
    entry = dict([('id', ''), ('header', ''), ('seq', '')])
    count = 0
    line_num = 0

    for line in fasta_file_f:
        line = line.strip()
        if line and line[0] == '>':
            count += 1
            key = '||'.join([entry[i] for i in keys if flags[i]])
            if key: # key != ''
                if key in fasta_dict: # check for duplicate key
                    logger.warning('%s : Line %d : Duplicate %s [%s] : ID = [%s].', fasta_file_f.name, line_num, '||'.join([i for i in keys if flags[i]]), key[:25] + (key[25:] and '..'), entry['id'])
                entry['seq'] = ''.join(entry['seq'])
                fasta_dict[key] = entry
                # check for url escaped id
                if id:
                    unescaped_id = unquote(entry['id'])
                    if id != unescaped_id:
                        key = '||'.join([unescaped_id] + [entry[i] for i in keys if i != 'id' and flags[i]])
                        entry['unescaped_id'] = unescaped_id
                        fasta_dict[key] = entry
                entry = dict()
            entry['header'] = line
            entry['id'] = line.split()[0][1:]
            entry['seq'] = []
        else:
            try:
                entry['seq'].append(line.upper())
            except:
                pass
        line_num += 1

    if isinstance(fasta_file, str):
        fasta_file_f.close()

    key = '||'.join([entry[i] for i in keys if flags[i]])
    if key: # key != ''
        if key in fasta_dict:
            logger.warning('%s : Line %d : Duplicate %s [%s] : ID = [%s].', fasta_file_f.name, line_num, '||'.join([i for i in keys if flags[i]]), key[:25] + (key[25:] and '..'), entry['id'])
        entry['seq'] = ''.join(entry['seq'])
        fasta_dict[key] = entry
        # check for url escaped id
        if id:
            unescaped_id = unquote(entry['id'])
            if id != unescaped_id:
                key = '||'.join([unescaped_id] + [entry[i] for i in keys if i != 'id' and flags[i]])
                entry['unescaped_id'] = unescaped_id
                fasta_dict[key] = entry

    return fasta_dict, count

def fasta_dict_to_file(fasta_dict, fasta_file, line_char_limit=None):
    """Write fasta_dict to fasta_file

    :param fasta_dict: returned by fasta_file_to_dict
    :param fasta_file: output file can be a string path or a file object
    :param line_char_limit: None = no limit (default)
    :return: None
    """
    fasta_fp = fasta_file
    if isinstance(fasta_file, str):
        fasta_fp = open(fasta_file, 'wb')

    for key in fasta_dict:
        seq = fasta_dict[key]['seq']
        if line_char_limit:
            seq = '\n'.join([seq[i:i+line_char_limit] for i in range(0, len(seq), line_char_limit)])
        fasta_fp.write(u'{0:s}\n{1:s}\n'.format(fasta_dict[key]['header'], seq))


class Gff3(object):
    def __init__(self, gff_file=None, fasta_external=None, logger=logger):
        self.logger = logger
        self.lines = []
        self.features = {}
        self.unresolved_parents = {}
        self.fasta_embedded = {}
        self.fasta_external = {}
        if gff_file:
            self.parse(gff_file)
        if fasta_external:
            self.parse_fasta_external(fasta_external)

    error_format = 'Line {current_line_num}: {error_type}: {message}\n-> {line}'

    def collect_descendants(self, line_data):
        collected_list = []
        if 'children' in line_data:
            children = line_data['children']
            for child in children:
                collected_list.append(child)
                collected_list.extend(self.collect_descendants(child))
        else:
            return
        return(collected_list)

    def collect_roots(self, line_data):
        collected_list = []
        try:
            if 'Parent' in line_data['attributes']:
                for p_line in line_data['parents']:
                    for p in p_line:
                        try:
                            if not 'Parent' in p['attributes']:
                                collected_list.append(p)
                            else:
                                collected_list.extend(self.collect_roots(p))
                        except:
                            pass
            else:
                collected_list.append(line_data)
                return(collected_list)
        except:
            pass
        return(collected_list)


    def add_line_error(self, line_data, error_info, log_level=logging.ERROR):
        """Helper function to record and log an error message

        :param line_data: dict
        :param error_info: dict
        :param logger:
        :param log_level: int
        :return:
        """
        if not error_info: return
        try:
            line_data['line_errors'].append(error_info)
        except KeyError:
            line_data['line_errors'] = [error_info]
        except TypeError: # no line_data
            pass
        try:
            self.logger.log(log_level, Gff3.error_format.format(current_line_num=line_data['line_index'] + 1, error_type=error_info['error_type'], message=error_info['message'], line=line_data['line_raw'].rstrip()))
        except AttributeError: # no logger
            pass

    def check_unresolved_parents(self):
        # check if any unresolved parents are now resolvable
        if len(self.unresolved_parents) > 0:
            self.logger.info('%d unresolved forward referencing parent ids, trying global lookup...' % len(self.unresolved_parents))
            globally_resolved_parents = set()
            for feature_id in self.unresolved_parents:
                if feature_id in self.features:
                    self.logger.info('  Resolved parent id: {0:s}, defined in lines: {1:s}, referenced in lines: {2:s}'.format(
                        feature_id,
                        ','.join([str(line_data['line_index'] + 1) for line_data in self.features[feature_id]]),
                        ','.join([str(line_data['line_index'] + 1) for line_data in self.unresolved_parents[feature_id]])))
                    globally_resolved_parents.add(feature_id)
                    for line_data in self.unresolved_parents[feature_id]:
                        line_data['parents'].append(self.features[feature_id])
                        for ld in self.features[feature_id]:
                            # no need to check if line_data in ld['children'], because it is impossible, each ld maps to only one feature_id, so the ld we get are all different
                            ld['children'].append(line_data)
            still_unresolved_parents = sorted(list(set(self.unresolved_parents) - globally_resolved_parents))
            if len(still_unresolved_parents) > 0:
                self.logger.info('{0:d} unresolved parent ids:'.format(len(still_unresolved_parents)))
                for feature_id in still_unresolved_parents:
                    self.logger.info('  Unresolved parent id: {0:s}, referenced in lines: {1:s}'.format(feature_id, ','.join(
                        [str(line_data['line_index'] + 1) for line_data in self.unresolved_parents[feature_id]])))

    def check_parent_boundary(self):
        """
        checks whether child features are within the coordinate boundaries of parent features

        :return:
        """
        check = 0
        flag = True
        for line in self.lines:
            if 'attributes' in line and 'ID' not in line['attributes'] and line['type'] in IDrequired:
                logger.error('[Missing ID] A model needs to have a unique ID, but this feature does not. Please fix it before running the program.\n\t\t- Line {0:s}: {1:s}'.format(str(line['line_index']+1), line['line_raw']))
                check += 1
        if check > 0:
            flag = False
            return flag

        for line in self.lines:
            for parent_feature in line['parents']:
                if len(parent_feature):
                    ok = False
                    for parent_line in parent_feature:
                        if parent_line['start'] <= line['start'] and line['end'] <= parent_line['end']:
                            ok = True
                            break
                    if not ok:
                        try:
                            self.add_line_error(line, {'message': '{2:s}: {0:s}: {1:s}'.format(parent_feature[0]['attributes']['ID'], ','.join(['({0:s}, {1:d}, {2:d})'.format(line['seqid'], line['start'], line['end']) for line in parent_feature]), ERROR_INFO['Ema0003']), 'error_type': 'BOUNDS', 'location': 'parent_boundary', 'eCode': 'Ema0003', 'error_level': 'Warning'}, log_level=logging.WARNING)
                        except:
                            logger.warning('Fail to check the boundary relationship between parent and child features (Ema0003)...\n\t- Line {0:s}: {1:s}'.format(str(line['line_index']+1), line['line_raw']))
                else:
                    logger.warning('[Parent feature missing] Cannot find the parents of this feature. Fail to check the boundary relationship (Ema0003).\n\t- Line {0:s}: {1:s}'.format(str(line['line_index']+1), line['line_raw']))
        return flag

    def check_phase(self, initial_phase):
        """
        1. get a list of CDS with the same parent
        2. sort according to strand
        3. calculate and validate phase
        """
        plus_minus = set(['+', '-'])
        for k, g in groupby(sorted([line for line in self.lines if  line['line_type'] == 'feature' and line['type'] == 'CDS' and 'Parent' in line['attributes']], key=lambda x: x['attributes']['Parent']), key=lambda x: x['attributes']['Parent']):
            cds_list = list(g)
            strand_set = list(set([line['strand'] for line in cds_list]))
            if len(strand_set) != 1:
                for line in cds_list:
                    self.add_line_error(line, {'message': 'Inconsistent CDS strand with parent: {0:s}'.format(str(k)), 'error_type': 'STRAND', 'eCode': 'Ema0007', 'error_level': 'Warning'}, log_level=logging.WARNING)
                continue
            if initial_phase:
                if len(cds_list) == 1:
                    if cds_list[0]['phase'] != 0:
                        self.add_line_error(cds_list[0], {'message': '{0:s} {1:d}, should be {2:d}'.format(ERROR_INFO['Ema0006'], cds_list[0]['phase'], 0), 'error_type': 'PHASE', 'eCode': 'Ema0006', 'error_level': 'Info'}, log_level=logging.INFO)
                    if type(cds_list[0]['phase']) != int:
                        logger.warning('[Phase] check_phase failed. \n\t\t- Line {0:s}: {1:s}'.format(str(cds_list[0]['line_index']+1), cds_list[0]['line_raw']))
                    continue
            strand = strand_set[0]
            if strand not in plus_minus:
                # don't process unknown strands
                continue
            if strand == '-':
                # sort end descending
                sorted_cds_list = sorted(cds_list, key=lambda x: x['end'], reverse=True)
            else:
                sorted_cds_list = sorted(cds_list, key=lambda x: x['start'])
            #If initial_phase is given, then program will test whether phase of the first CDS is 0. If it is not, then show warning.
            #If initial_phase is not given (default), then program will not test whether phase of the first CDS is 0
            if initial_phase:
                phase = 0
            else:
                phase = sorted_cds_list[0]['phase']
                if type(phase) != int:
                    logger.warning('[Phase] check_phase failed. \n\t\t- Line {0:s}: {1:s}'.format(str(sorted_cds_list[0]['line_index']+1), sorted_cds_list[0]['line_raw']))
            for line in sorted_cds_list:
                if line['phase'] != phase:
                    try:
                        self.add_line_error(line, {'message': 'Wrong phase {0:d}, should be {1:d}'.format(line['phase'], phase), 'error_type': 'PHASE', 'eCode': 'Ema0006', 'error_level': 'Info'}, log_level=logging.INFO)
                    except:
                        logger.warning('[Phase] check_phase failed. \n\t\t- Line {0:s}: {1:s}'.format(str(line['line_index']+1), line['line_raw']))
                try:
                    phase = (3 - ((line['end'] - line['start'] + 1 - phase) % 3)) % 3
                except:
                    pass
    def parse_fasta_external(self, fasta_file):
        self.fasta_external, count = fasta_file_to_dict(fasta_file)

    def check_reference(self, sequence_region=False, fasta_embedded=False, fasta_external=False, check_bounds=True, check_n=True, allowed_num_of_n=0, feature_types=('CDS',)):
        """
        Check seqid, bounds and the number of Ns in each feature using one or more reference sources.

        Seqid check: check if the seqid can be found in the reference sources.

        Bounds check: check the start and end fields of each features and log error if the values aren't within the seqid sequence length, requires at least one of these sources: ##sequence-region, embedded #FASTA, or external FASTA file.

        Ns check: count the number of Ns in each feature with the type specified in *line_types (default: 'CDS') and log an error if the number is greater than allowed_num_of_n (default: 0), requires at least one of these sources: embedded #FASTA, or external FASTA file.

        When called with all source parameters set as False (default), check all available sources, and log debug message if unable to perform a check due to none of the reference sources being available.

        If any source parameter is set to True, check only those sources marked as True, log error if those sources don't exist.

        :param sequence_region: check bounds using the ##sequence-region directive (default: False)
        :param fasta_embedded: check bounds using the embedded fasta specified by the ##FASTA directive (default: False)
        :param fasta_external: check bounds using the external fasta given by the self.parse_fasta_external (default: False)
        :param check_bounds: If False, don't run the bounds check (default: True)
        :param check_n: If False, don't run the Ns check (default: True)
        :param allowed_num_of_n: only report features with a number of Ns greater than the specified value (default: 0)
        :param feature_types: only check features of these feature_types, multiple types may be specified, if none are specified, check only 'CDS'
        :return: error_lines: a set of line_index(int) with errors detected by check_reference
        """
        # collect lines with errors in this set
        error_lines = set()
        # check if we have a parsed gff3
        if not self.lines:
            self.logger.debug('.parse(gff_file) before calling .check_bounds()')
            return error_lines
        # setup default line_types
        check_n_feature_types = set(feature_types)
        if len(check_n_feature_types) == 0:
            check_n_feature_types.add('CDS')
        # compile regex
        n_segments_finditer = re.compile(r'[Nn]+').finditer
        # check_all_sources mode
        check_all_sources = True
        if sequence_region or fasta_embedded or fasta_external:
            check_all_sources = False
        # get a list of line_data with valid start and end coordinates and unescape the seqid
        start_end_error_locations = set(('start', 'end', 'start,end'))
        valid_line_data_seqid = [(line_data, unquote(line_data['seqid'])) for line_data in self.lines if line_data['line_type'] == 'feature' and line_data['seqid'] != '.' and (not line_data['line_errors'] or not [error_info for error_info in line_data['line_errors'] if 'location' in error_info and error_info['location'] in start_end_error_locations])]
        checked_at_least_one_source = False
        # check directive
        # don't use any directives with errors
        valid_sequence_regions = dict([(unquote(line_data['seqid']), line_data) for line_data in self.lines if line_data['directive'] == '##sequence-region' and not line_data['line_errors']])
        unresolved_seqid = set()
        if (check_all_sources or sequence_region) and valid_sequence_regions:
            checked_at_least_one_source = True
            for line_data, seqid in valid_line_data_seqid:
                if seqid not in valid_sequence_regions and seqid not in unresolved_seqid:
                    unresolved_seqid.add(seqid)
                    error_lines.add(line_data['line_index'])
                    self.add_line_error(line_data, {'message': 'Seqid not found in any ##sequence-region: {0:s}'.format(
                        seqid), 'error_type': 'BOUNDS', 'location': 'sequence_region', 'eCode': 'Esf0004'})
                    continue
                if seqid not in unresolved_seqid:
                    if line_data['start'] < valid_sequence_regions[seqid]['start']:
                        error_lines.add(line_data['line_index'])
                        self.add_line_error(line_data, {'message': 'Start is less than the ##sequence-region start: %d' % valid_sequence_regions[seqid]['start'], 'error_type': 'BOUNDS', 'location': 'sequence_region', 'eCode': 'Esf0005'})
                    if line_data['end'] > valid_sequence_regions[seqid]['end']:
                        error_lines.add(line_data['line_index'])
                        self.add_line_error(line_data, {'message': 'End is greater than the ##sequence-region end: %d' % valid_sequence_regions[seqid]['end'], 'error_type': 'BOUNDS', 'location': 'sequence_region', 'eCode': 'Esf0006'})
        elif sequence_region:
            self.logger.debug('##sequence-region not found in GFF3')
        # check fasta_embedded
        unresolved_seqid = set()
        if (check_all_sources or fasta_embedded) and self.fasta_embedded:
            checked_at_least_one_source = True
            for line_data, seqid in valid_line_data_seqid:
                if seqid not in self.fasta_embedded and seqid not in unresolved_seqid:
                    unresolved_seqid.add(seqid)
                    error_lines.add(line_data['line_index'])
                    self.add_line_error(line_data, {'message': 'Seqid not found in the embedded ##FASTA: %s' % seqid, 'error_type': 'BOUNDS', 'location': 'fasta_embedded', 'eCode': 'Esf0007'})
                    continue
                # check bounds
                if seqid not in unresolved_seqid:
                    if line_data['end'] > len(self.fasta_embedded[seqid]['seq']):
                        error_lines.add(line_data['line_index'])
                        self.add_line_error(line_data, {'message': 'End is greater than the embedded ##FASTA sequence length: %d' % len(self.fasta_embedded[seqid]['seq']), 'error_type': 'BOUNDS', 'location': 'fasta_embedded', 'eCode': 'Esf0008'})
                # check n
                if check_n and line_data['type'] in check_n_feature_types:
                    """
                    >>> timeit("a.lower().count('n')", "import re; a = ('ASDKADSJHFIUDNNNNNNNnnnnSHFD'*50)")
                    5.540903252684302
                    >>> timeit("a.count('n'); a.count('N')", "import re; a = ('ASDKADSJHFIUDNNNNNNNnnnnSHFD'*50)")
                    2.3504867946058425
                    >>> timeit("re.findall('[Nn]+', a)", "import re; a = ('ASDKADSJHFIUDNNNNNNNnnnnSHFD'*50)")
                    30.60731204915959
                    """
                    n_count = self.fasta_embedded[seqid]['seq'].count('N', line_data['start'] - 1, line_data['end']) + self.fasta_embedded[seqid]['seq'].count('n', line_data['start'] - 1, line_data['end'])
                    if n_count > allowed_num_of_n:
                        # get detailed segments info
                        n_segments = [(m.start(), m.end() - m.start()) for m in n_segments_finditer(self.fasta_embedded[seqid]['seq'], line_data['start'] - 1, line_data['end'])]
                        n_segments_str = ['(%d, %d)' % (m[0], m[1]) for m in n_segments]
                        error_lines.add(line_data['line_index'])
                        self.add_line_error(line_data, {'message': 'Found %d Ns in %s feature of length %d using the embedded ##FASTA, consists of %d segment (start, length): %s' % (n_count, line_data['type'], line_data['end'] - line_data['start'], len(n_segments), ', '.join(n_segments_str)), 'error_type': 'N_COUNT', 'n_segments': n_segments, 'location': 'fasta_embedded', 'eCode': 'Esf0009', 'error_level': 'Info'}, log_level=logging.INFO)

        elif fasta_embedded:
            self.logger.debug('Embedded ##FASTA not found in GFF3')
        # check fasta_external
        unresolved_seqid = set()
        if (check_all_sources or fasta_external) and self.fasta_external:
            checked_at_least_one_source = True
            for line_data, seqid in valid_line_data_seqid:
                if seqid not in self.fasta_external and seqid not in unresolved_seqid:
                    unresolved_seqid.add(seqid)
                    error_lines.add(line_data['line_index'])
                    self.add_line_error(line_data, {'message': 'Seqid not found in the external FASTA file: %s' % seqid, 'error_type': 'BOUNDS', 'location': 'fasta_external', 'eCode': 'Esf0010'})
                    continue
                # check bounds
                if seqid not in unresolved_seqid:
                    if line_data['end'] > len(self.fasta_external[seqid]['seq']):
                        error_lines.add(line_data['line_index'])
                        self.add_line_error(line_data, {'message': 'End is greater than the external FASTA sequence length: %d' % len(self.fasta_external[seqid]['seq']), 'error_type': 'BOUNDS', 'location': 'fasta_external', 'eCode': 'Esf0011'})
                # check n
                if check_n and line_data['type'] in check_n_feature_types:
                    try:
                        n_count = self.fasta_external[seqid]['seq'].count('N', line_data['start'] - 1, line_data['end']) + self.fasta_external[seqid]['seq'].count('n', line_data['start'] - 1, line_data['end'])
                        if n_count > allowed_num_of_n:
                            # get detailed segments info
                            n_segments = [(m.start(), m.end() - m.start()) for m in n_segments_finditer(self.fasta_external[seqid]['seq'], line_data['start'] - 1, line_data['end'])]
                            n_segments_str = ['(%d, %d)' % (m[0], m[1]) for m in n_segments]
                            error_lines.add(line_data['line_index'])
                            self.add_line_error(line_data, {'message': 'Found %d Ns in %s feature of length %d using the external FASTA, consists of %d segment (start, length): %s' % (n_count, line_data['type'], line_data['end'] - line_data['start'], len(n_segments), ', '.join(n_segments_str)), 'error_type': 'N_COUNT', 'n_segments': n_segments, 'location': 'fasta_external', 'eCode': 'Esf0012', 'error_level': 'Info'}, log_level=logging.INFO)
                    except:
                        logger.warning('Sequence ID {0:s} not found in FASTA file.'.format(seqid))
        elif fasta_external:
            self.logger.debug('External FASTA file not given')
        if check_all_sources and not checked_at_least_one_source:
            self.logger.debug('Unable to perform bounds check, requires at least one of the following sources: ##sequence-region, embedded ##FASTA, or external FASTA file')
        return error_lines

    def parse(self, gff_file, strict=False):
        """Parse the gff file into the following data structures:

        * lines(list of line_data(dict))
            - line_index(int): the index in lines
            - line_raw(str)
            - line_type(str in ['feature', 'directive', 'comment', 'blank', 'unknown'])
            - line_errors(list of str): a list of error messages
            - line_status(str in ['normal', 'modified', 'removed'])
            - parents(list of feature(list of line_data(dict))): may have multiple parents
            - children(list of line_data(dict))
            - extra fields depending on line_type
            * directive
                - directive(str in ['##gff-version', '##sequence-region', '##feature-ontology', '##attribute-ontology', '##source-ontology', '##species', '##genome-build', '###', '##FASTA'])
                - extra fields depending on directive
            * feature
                - seqid(str): must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|] using RFC 3986 Percent-Encoding
                - source(str)
                - type(str in so_types)
                - start(int)
                - end(int)
                - score(float)
                - strand(str in ['+', '-', '.', '?'])
                - phase(int in [0, 1, 2])
                - attributes(dict of tag(str) to value)
                    - ID(str)
                    - Name(str)
                    - Alias(list of str): multi value
                    - Parent(list of str): multi value
                    - Target(dict)
                        - target_id(str)
                        - start(int)
                        - end(int)
                        - strand(str in ['+', '-', ''])
                    - Gap(str): CIGAR format
                    - Derives_from(str)
                    - Note(list of str): multi value
                    - Dbxref(list of str): multi value
                    - Ontology_term(list of str): multi value
                    - Is_circular(str in ['true'])
            * fasta_dict(dict of id(str) to sequence_item(dict))
                - id(str)
                - header(str)
                - seq(str)
                - line_length(int)

        * features(dict of feature_id(str in line_data['attributes']['ID']) to feature(list of line_data(dict)))

        A feature is a list of line_data(dict), since all lines that share an ID collectively represent a single feature.

        During serialization, line_data(dict) references should be converted into line_index(int)

        :param gff_file: a string path or file object
        :param strict: when true, throw exception on syntax and format errors. when false, use best effort to finish parsing while logging errors
        """
        valid_strand = set(('+', '-', '.', '?'))
        valid_phase = set((0, 1, 2))
        multi_value_attributes = set(('replace','Parent', 'Alias', 'Note', 'Dbxref', 'Ontology_term'))
        valid_attribute_target_strand = set(('+', '-', ''))
        reserved_attributes = set(('replace','ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from', 'Note', 'Dbxref', 'Ontology_term', 'Is_circular'))

        # illegal character check
        # Literal use of tab, newline, carriage return, the percent (%) sign, and control characters must be encoded using RFC 3986 Percent-Encoding; no other characters may be encoded.
        # control characters: \x00-\x1f\x7f this includes tab(\x09), newline(\x0a), carriage return(\x0d)
        # seqid may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]
        # URL escaping rules are used for tags or values containing the following characters: ",=;".
        #>>> timeit("unescaped_seqid('Un.7589')", "import re; unescaped_seqid = re.compile(r'[^a-zA-Z0-9.:^*$@!+_?|%-]|%(?![0-9a-fA-F]{2})').search")
        #0.4128372745785036
        #>>> timeit("unescaped_seqid2('Un.7589')", "import re; unescaped_seqid2 = re.compile(r'^([a-zA-Z0-9.:^*$@!+_?|-]|%[0-9a-fA-F]{2})+$').search")
        #0.9012313532265175
        unescaped_seqid = re.compile(r'[^a-zA-Z0-9.:^*$@!+_?|%-]|%(?![0-9a-fA-F]{2})').search
        unescaped_field = re.compile(r'[\x00-\x1f\x7f]|%(?![0-9a-fA-F]{2})').search

        gff_fp = gff_file
        if isinstance(gff_file, str):
            if gff_file.endswith('.gz'):
                gff_fp = gzip.open(gff_file, 'rt')  # gzip support
            else:
                gff_fp = open(gff_file, 'r')

        lines = []
        current_line_num = 1 # line numbers start at 1
        features = defaultdict(list)
        # key = the unresolved id, value = a list of line_data(dict)
        unresolved_parents = defaultdict(list)

        for line_raw in gff_fp:
            line_data = {
                'line_index': current_line_num - 1,
                'line_raw': line_raw,
                'line_status': 'normal',
                'parents': [],
                'children': [],
                'line_type': '',
                'directive': '',
                'line_errors': [],
                'type': '',
            }
            line_strip = line_raw.strip()
            if line_strip != line_raw[:len(line_strip)]:
                self.add_line_error(line_data, {'message': 'White chars not allowed at the start of a line', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0013'})
            if current_line_num == 1 and not line_strip.startswith('##gff-version'):
                self.add_line_error(line_data, {'message': '"##gff-version" missing from the first line', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0014'})
            if len(line_strip) == 0:
                line_data['line_type'] = 'blank'
                current_line_num += 1
                lines.append(line_data)
                continue
            if line_strip.startswith('##'):
                line_data['line_type'] = 'directive'
                if line_strip.startswith('##sequence-region'):
                    # ##sequence-region seqid start end
                    # This element is optional, but strongly encouraged because it allows parsers to perform bounds checking on features.
                    # only one ##sequence-region directive may be given for any given seqid
                    # all features on that landmark feature (having that seqid) must be contained within the range defined by that ##sequence-region diretive. An exception to this rule is allowed when a landmark feature is marked with the Is_circular attribute.
                    line_data['directive'] = '##sequence-region'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 3:
                        self.add_line_error(line_data, {'message': 'Expecting 3 fields, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['seqid'] = tokens[0]
                        # check for duplicate ##sequence-region seqid
                        if [True for d in lines if ('directive' in d and d['directive'] == '##sequence-region' and 'seqid' in d and d['seqid'] == line_data['seqid'])]:
                            self.add_line_error(line_data, {'message': '##sequence-region seqid: "%s" may only appear once' % line_data['seqid'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0016'})
                        try:
                            all_good = True
                            try:
                                line_data['start'] = int(tokens[1])
                                if line_data['start'] < 1:
                                    self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0002'], tokens[1]), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0002'})
                            except ValueError:
                                all_good = False
                                self.add_line_error(line_data, {'message': 'Start is not a valid integer: "%s"' % tokens[1], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0017'})
                                line_data['start'] = tokens[1]
                            try:
                                line_data['end'] = int(tokens[2])
                                if line_data['end'] < 1:
                                    self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0002'], tokens[2]), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0002'})
                            except ValueError:
                                all_good = False
                                self.add_line_error(line_data, {'message': 'End is not a valid integer: "%s"' % tokens[2], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0017'})
                                line_data['start'] = tokens[2]
                            # if all_good then both start and end are int, so we can check if start is not less than or equal to end
                            if all_good and line_data['start'] > line_data['end']:
                                self.add_line_error(line_data, {'message': 'Start is not less than or equal to end', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0018'})
                        except IndexError:
                            pass
                elif line_strip.startswith('##gff-version'):
                    # The GFF version, always 3 in this specification must be present, must be the topmost line of the file and may only appear once in the file.
                    line_data['directive'] = '##gff-version'
                    # check if it appeared before
                    if [True for d in lines if ('directive' in d and d['directive'] == '##gff-version')]:
                        self.add_line_error(line_data, {'message': '##gff-version missing from the first line', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0014'})
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 1:
                        self.add_line_error(line_data, {'message': 'Expecting 1 field, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        try:
                            line_data['version'] = int(tokens[0])
                            if line_data['version'] != 3:
                                self.add_line_error(line_data, {'message': 'Version is not "3": "%s"' % tokens[0], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0019', 'error_level': 'Info'}, log_level=logging.INFO)
                        except ValueError:
                            self.add_line_error(line_data, {'message': 'Version is not a valid integer: "%s"' % tokens[0], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0020'})
                            line_data['version'] = tokens[0]
                elif line_strip.startswith('###'):
                    # This directive (three # signs in a row) indicates that all forward references to feature IDs that have been seen to this point have been resolved.
                    line_data['directive'] = '###'
                elif line_strip.startswith('##FASTA'):
                    # This notation indicates that the annotation portion of the file is at an end and that the
                    # remainder of the file contains one or more sequences (nucleotide or protein) in FASTA format.
                    line_data['directive'] = '##FASTA'
                    self.logger.info('Reading embedded ##FASTA sequence')
                    self.fasta_embedded, count = fasta_file_to_dict(gff_fp)
                    self.logger.info('%d sequences read' % len(self.fasta_embedded))
                elif line_strip.startswith('##feature-ontology'):
                    # ##feature-ontology URI
                    # This directive indicates that the GFF3 file uses the ontology of feature types located at the indicated URI or URL.
                    line_data['directive'] = '##feature-ontology'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 1:
                        self.add_line_error(line_data, {'message': 'Expecting 1 field, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['URI'] = tokens[0]
                elif line_strip.startswith('##attribute-ontology'):
                    # ##attribute-ontology URI
                    # This directive indicates that the GFF3 uses the ontology of attribute names located at the indicated URI or URL.
                    line_data['directive'] = '##attribute-ontology'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 1:
                        self.add_line_error(line_data, {'message': 'Expecting 1 field, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['URI'] = tokens[0]
                elif line_strip.startswith('##source-ontology'):
                    # ##source-ontology URI
                    # This directive indicates that the GFF3 uses the ontology of source names located at the indicated URI or URL.
                    line_data['directive'] = '##source-ontology'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 1:
                        self.add_line_error(line_data, {'message': 'Expecting 1 field, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['URI'] = tokens[0]
                elif line_strip.startswith('##species'):
                    # ##species NCBI_Taxonomy_URI
                    # This directive indicates the species that the annotations apply to.
                    line_data['directive'] = '##species'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 1:
                        self.add_line_error(line_data, {'message': 'Expecting 1 field, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['NCBI_Taxonomy_URI'] = tokens[0]
                elif line_strip.startswith('##genome-build'):
                    # ##genome-build source buildName
                    # The genome assembly build name used for the coordinates given in the file.
                    line_data['directive'] = '##genome-build'
                    tokens = line_strip.split()[1:]
                    if len(tokens) != 2:
                        self.add_line_error(line_data, {'message': 'Expecting 2 fields, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0015'})
                    if len(tokens) > 0:
                        line_data['source'] = tokens[0]
                        try:
                            line_data['buildName'] = tokens[1]
                        except IndexError:
                            pass
                else:
                    self.add_line_error(line_data, {'message': 'Unknown directive', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0021', 'error_level': 'Info'}, log_level=logging.INFO)
                    tokens = line_strip.split()
                    line_data['directive'] = tokens[0]
            elif line_strip.startswith('#'):
                line_data['line_type'] = 'comment'
            else:
                # line_type may be a feature or unknown
                line_data['line_type'] = 'feature'
                tokens =  list(map(str.strip, line_raw.split('\t')))
                if len(tokens) != 9:
                    self.add_line_error(line_data, {'message': 'Features should contain 9 fields, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0022'})
                for i, t in enumerate(tokens):
                    if not t:
                        self.add_line_error(line_data, {'message': 'Empty field: %d, must have a "."' % (i + 1), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0022'})
                try:
                    line_data['seqid'] = tokens[0]
                    if unescaped_seqid(tokens[0]):
                        self.add_line_error(line_data, {'message': 'Seqid must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]: "%s"' % tokens[0], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0023'})
                    line_data['source'] = tokens[1]
                    if unescaped_field(tokens[1]):
                        self.add_line_error(line_data, {'message': 'Source must escape the percent (%%) sign and any control characters: "%s"' % tokens[1], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0023'})
                    line_data['type'] = tokens[2]
                    if unescaped_field(tokens[2]):
                        self.add_line_error(line_data, {'message': 'Type must escape the percent (%%) sign and any control characters: "%s"' % tokens[2], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0023'})
                    all_good = True
                    try:
                        line_data['start'] = int(tokens[3])
                        if line_data['start'] < 1:
                            self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0002'], tokens[3]), 'error_type': 'FORMAT', 'location': 'start', 'eCode': 'Esf0002'})
                    except ValueError:
                        all_good = False
                        line_data['start'] = tokens[3]
                        if line_data['start'] != '.':
                            self.add_line_error(line_data, {'message': 'Start is not a valid integer: "%s"' % line_data['start'], 'error_type': 'FORMAT', 'location': 'start', 'eCode': 'Esf0017'})
                    try:
                        line_data['end'] = int(tokens[4])
                        if line_data['end'] < 1:
                            self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0002'], tokens[4]), 'error_type': 'FORMAT', 'location': 'end', 'eCode': 'Esf0002'})
                    except ValueError:
                        all_good = False
                        line_data['end'] = tokens[4]
                        if line_data['end'] != '.':
                            self.add_line_error(line_data, {'message': 'End is not a valid integer: "%s"' % line_data['end'], 'error_type': 'FORMAT', 'location': 'end', 'eCode': 'Esf0017'})
                    # if all_good then both start and end are int, so we can check if start is not less than or equal to end
                    if all_good and line_data['start'] > line_data['end']:
                        self.add_line_error(line_data, {'message': 'Start is not less than or equal to end', 'error_type': 'FORMAT', 'location': 'start,end', 'eCode': 'Esf0018'})
                    try:
                        line_data['score'] = float(tokens[5])
                    except ValueError:
                        line_data['score'] = tokens[5]
                        if line_data['score'] != '.':
                            self.add_line_error(line_data, {'message': 'Score is not a valid floating point number: "%s"' % line_data['score'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0024'})
                    line_data['strand'] = tokens[6]
                    if line_data['strand'] not in valid_strand: # set(['+', '-', '.', '?'])
                        self.add_line_error(line_data, {'message': 'Strand has illegal characters: "%s"' % tokens[6], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0025'})
                    try:
                        line_data['phase'] = int(tokens[7])
                        if line_data['phase'] not in valid_phase: # set([0, 1, 2])
                            self.add_line_error(line_data, {'message': 'Phase is not 0, 1, or 2: "%s"' % tokens[7], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0026'})
                    except ValueError:
                        line_data['phase'] = tokens[7]
                        if line_data['phase'] != '.':
                            self.add_line_error(line_data, {'message': 'Phase is not a valid integer: "%s"' % line_data['phase'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0026'})
                        elif line_data['type'] == 'CDS':
                            self.add_line_error(line_data, {'message': 'Phase is required for all CDS features', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0027'})
                    # parse attributes, ex: ID=exon00003;Parent=mRNA00001,mRNA00003;Name=EXON.1
                    # URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape.
                    # Note that attribute names are case sensitive. "Parent" is not the same as "parent".
                    # All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with a lowercase letter can be used freely by applications.
                    if unescaped_field(tokens[8]):
                        self.add_line_error(line_data, {'message': 'Attributes must escape the percent (%) sign and any control characters', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0028', 'error_level': 'Info'}, log_level=logging.INFO)
                    attribute_tokens = tuple(tuple(t for t in a.split('=')) for a in tokens[8].split(';') if a)
                    line_data['attributes'] = {}
                    if len(attribute_tokens) == 1 and len(attribute_tokens[0]) == 1 and attribute_tokens[0][0] == '.':
                        pass # no attributes
                    else:
                        for a in attribute_tokens:
                            if len(a) != 2:
                                self.add_line_error(line_data, {'message': 'Attributes must contain one and only one equal (=) sign: "%s"' % ('='.join(a)), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0029'})
                            try:
                                tag, value = a
                            except ValueError:
                                tag, value = a[0], ''
                            if not tag:
                                self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0030'], '='.join(a)), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0030'})
                            if not value.strip():
                                self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0031'], '='.join(a)), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0031'})
                            if tag in line_data['attributes']:
                                self.add_line_error(line_data, {'message': '%s: "%s"' % (ERROR_INFO['Esf0032'], tag), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0032'})
                            if tag in multi_value_attributes: # set(['replace', 'Parent', 'Alias', 'Note', 'Dbxref', 'Ontology_term'])
                                if value.find(', ') >= 0 or value.find(' ,') >= 0:
                                    self.add_line_error(line_data, {'message': 'Found ", " in %s attribute, possible unescaped ",": "%s"' % (tag, value), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0033', 'error_level': 'Info'}, log_level=logging.INFO)
                                # In addition to Parent, the Alias, Note, Dbxref and Ontology_term attributes can have multiple values.
                                if tag in line_data['attributes']: # if this tag has been seen before
                                    if tag == 'Note': # don't check for duplicate notes
                                        line_data['attributes'][tag].extend(value.split(','))
                                    else: # only add non duplicate values
                                        line_data['attributes'][tag].extend([s for s in value.split(',') if s not in line_data['attributes'][tag]])
                                else:
                                    line_data['attributes'][tag] = value.split(',')
                                # check for duplicate values
                                if tag != 'Note' and len(line_data['attributes'][tag]) != len(set(line_data['attributes'][tag])):
                                    count_values = [(len(list(group)), key) for key, group in groupby(sorted(line_data['attributes'][tag]))]
                                    self.add_line_error(line_data, {'message': '%s %s: %s' % (tag, ERROR_INFO['Esf0034'],', '.join(['(%d, %s)' % (c, v) for c, v in count_values if c > 1])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0034', 'error_level': 'Info'}, log_level=logging.INFO)
                                    # remove duplicate
                                    line_data['attributes'][tag] = list(set(line_data['attributes'][tag]))

                                if tag == 'Parent':
                                    for feature_id in line_data['attributes']['Parent']:
                                        try:
                                            line_data['parents'].append(features[feature_id])
                                            for ld in features[feature_id]:
                                                # no need to check if line_data in ld['children'], because it is impossible, each ld maps to only one feature_id, so the ld we get are all different
                                                ld['children'].append(line_data)
                                        except KeyError: # features[id]
                                            self.add_line_error(line_data, {'message': '%s attribute has unresolved forward reference: %s' % (tag, feature_id), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0035', 'error_level': 'Info'}, log_level=logging.INFO)
                                            unresolved_parents[feature_id].append(line_data)
                            elif tag == 'Target':
                                if value.find(',') >= 0:
                                    self.add_line_error(line_data, {'message': 'Value of %s attribute contains unescaped ",": "%s"' % (tag, value), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0036', 'error_level': 'Info'}, log_level=logging.INFO)
                                target_tokens = value.split(' ')
                                if len(target_tokens) < 3 or len(target_tokens) > 4:
                                    self.add_line_error(line_data, {'message': 'Target attribute should have 3 or 4 values, got %d: %s' % (len(target_tokens), repr(tokens)), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0037'})
                                line_data['attributes'][tag] = {}
                                try:
                                    line_data['attributes'][tag]['target_id'] = target_tokens[0]
                                    all_good = True
                                    try:
                                        line_data['attributes'][tag]['start'] = int(target_tokens[1])
                                        if line_data['attributes'][tag]['start'] < 1:
                                            self.add_line_error(line_data, {'message': 'Start value of Target attribute is not a valid 1-based integer coordinate: "%s"' % target_tokens[1], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0038'})
                                    except ValueError:
                                        all_good = False
                                        line_data['attributes'][tag]['start'] = target_tokens[1]
                                        self.add_line_error(line_data, {'message': 'Start value of Target attribute is not a valid integer: "%s"' % line_data['attributes'][tag]['start'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0038'})
                                    try:
                                        line_data['attributes'][tag]['end'] = int(target_tokens[2])
                                        if line_data['attributes'][tag]['end'] < 1:
                                            self.add_line_error(line_data, {'message': 'End value of Target attribute is not a valid 1-based integer coordinate: "%s"' % target_tokens[2], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0038'})
                                    except ValueError:
                                        all_good = False
                                        line_data['attributes'][tag]['end'] = target_tokens[2]
                                        self.add_line_error(line_data, {'message': 'End value of Target attribute is not a valid integer: "%s"' % line_data['attributes'][tag]['end'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0038'})
                                    # if all_good then both start and end are int, so we can check if start is not less than or equal to end
                                    if all_good and line_data['attributes'][tag]['start'] > line_data['attributes'][tag]['end']:
                                        self.add_line_error(line_data, {'message': 'Start is not less than or equal to end', 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0018'})
                                    line_data['attributes'][tag]['strand'] = target_tokens[3]
                                    if line_data['attributes'][tag]['strand'] not in valid_attribute_target_strand: # set(['+', '-', ''])
                                        self.add_line_error(line_data, {'message': 'Strand value of Target attribute has illegal characters: "%s"' % line_data['attributes'][tag]['strand'], 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0039'})
                                except IndexError:
                                    pass
                            else:
                                if value.find(',') >= 0:
                                    self.add_line_error(line_data, {'message': 'Value of %s attribute contains unescaped ",": "%s"' % (tag, value), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0036'}, log_level=logging.INFO)
                                line_data['attributes'][tag] = value
                                if tag == 'Is_circular' and value != 'true':
                                    self.add_line_error(line_data, {'message': 'Value of Is_circular attribute is not "true": "%s"' % value, 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0040'})
                                elif tag[:1].isupper() and tag not in reserved_attributes: # {'replace','ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from', 'Note', 'Dbxref', 'Ontology_term', 'Is_circular'}
                                    self.add_line_error(line_data, {'message': 'Unknown reserved (uppercase) attribute: "%s"' % tag, 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0041'})
                                elif tag == 'ID':
                                    # check for duplicate ID in non-adjacent lines
                                    try:
                                        if value in features and 'attributes' in lines[-1] and lines[-1]['attributes'][tag] != value:
                                            self.add_line_error(line_data, {'message': '%s: "%s" in non-adjacent lines: %s' % (ERROR_INFO['Emr0003'], value, ','.join([str(f['line_index'] + 1) for f in features[value]])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Emr0003', 'error_level': 'Error'})
                                        elif value in features and 'attributes' not in lines[-1]:
                                            self.add_line_error(line_data, {'message': '%s: "%s" in non-adjacent lines: %s' % (ERROR_INFO['Emr0003'], value, ','.join([str(f['line_index'] + 1) for f in features[value]])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Emr0003', 'error_level': 'Error'})
                                    except:
                                        logger.warning('[Missing ID] Program failed. \n\t\t- Line {0:s}: {1:s}'.format(str(lines[-1]['line_index']+1), lines[-1]['line_raw']))

                                    features[value].append(line_data)
                except IndexError:
                    pass
            current_line_num += 1
            lines.append(line_data)

        if isinstance(gff_file, str):
            gff_fp.close()

        # global look up of unresolved parents
        for feature_id in unresolved_parents:
            if feature_id in features:
                for line in unresolved_parents[feature_id]:
                    self.add_line_error(line, {'message': 'Unresolved forward reference: "%s", found defined in lines: %s' % (feature_id, ','.join([str(ld['line_index'] + 1) for ld in features[feature_id]])), 'error_type': 'FORMAT', 'location': '', 'eCode': 'Esf0042'})

        self.lines = lines
        self.features = features
        return 1

    def descendants(self, line_data):
        """
        BFS graph algorithm
        :param line_data: line_data(dict) with line_data['line_index'] or line_index(int)
        :return: list of line_data(dict)
        """
        # get start node
        try:
            start = line_data['line_index']
        except TypeError:
            start = self.lines[line_data]['line_index']
        visited_set, visited_list, queue = set(), [], [start]
        while queue:
            node = queue.pop(0)
            if node not in visited_set:
                visited_set.add(node)
                visited_list.append(self.lines[node])
                #queue.extend([ld['line_index'] for ld in self.lines[node]['children'] if ld['line_index'] not in visited_set])
                ### To write out gff file follwing the order of gene, mRNA, exon, CDS
                #'''
                for ld in self.lines[node]['children']:
                    if ld['line_index'] not in visited_set:
                        queue.extend([ld['line_index']])
                    for gld in ld['children']:
                        if gld['line_index'] not in visited_set:
                            queue.extend([gld['line_index']])
                #'''
                ### To write out gff file follwing the order of gene, mRNA, exon, CDS
        return visited_list[1:]

    def ancestors(self, line_data):
        """
        BFS graph algorithm

        :param line_data: line_data(dict) with line_data['line_index'] or line_index(int)
        :return: list of line_data(dict)
        """
        # get start node
        try:
            start = line_data['line_index']
        except TypeError:
            start = self.lines[line_data]['line_index']
        visited_set, visited_list, queue = set(), [], [start]
        while queue:
            node = queue.pop(0)
            if node not in visited_set:
                visited_set.add(node)
                visited_list.append(self.lines[node])
                queue.extend([ld['line_index'] for f in self.lines[node]['parents'] for ld in f if ld['line_index'] not in visited_set])
        return visited_list[1:]

    def adopt(self, old_parent, new_parent):
        """
        Transfer children from old_parent to new_parent

        :param old_parent: feature_id(str) or line_index(int) or line_data(dict) or feature
        :param new_parent: feature_id(str) or line_index(int) or line_data(dict)
        :return: List of children transferred
        """
        try: # assume line_data(dict)
            old_id = old_parent['attributes']['ID']
        except TypeError:
            try: # assume line_index(int)
                old_id = self.lines[old_parent]['attributes']['ID']
            except TypeError: # assume feature_id(str)
                old_id = old_parent
        old_feature = self.features[old_id]
        old_indexes = [ld['line_index'] for ld in old_feature]
        try: # assume line_data(dict)
            new_id = new_parent['attributes']['ID']
        except TypeError:
            try: # assume line_index(int)
                new_id = self.lines[new_parent]['attributes']['ID']
            except TypeError: # assume feature_id(str)
                new_id = new_parent
        new_feature = self.features[new_id]
        new_indexes = [ld['line_index'] for ld in new_feature]
        # build a list of children to be moved
        # add the child to the new parent's children list if its not already there
        # update the child's parent list and parent attribute
        # finally remove the old parent's children list
        children = old_feature[0]['children']
        new_parent_children_set = set([ld['line_index'] for ld in new_feature[0]['children']])
        for child in children:
            if child['line_index'] not in new_parent_children_set:
                new_parent_children_set.add(child['line_index'])
                for new_ld in new_feature:
                    new_ld['children'].append(child)
                child['parents'].append(new_feature)
                child['attributes']['Parent'].append(new_id)
            # remove multiple, list.remove() only removes 1
            child['parents'] = [f for f in child['parents'] if f[0]['attributes']['ID'] != old_id]
            child['attributes']['Parent'] = [d for d in child['attributes']['Parent'] if d != old_id]
        for old_ld in old_feature:
            old_ld['children'] = []
        return children

    def adopted(self, old_child, new_child):
        """
        Transfer parents from old_child to new_child

        :param old_child: line_data(dict) with line_data['line_index'] or line_index(int)
        :param new_child: line_data(dict) with line_data['line_index'] or line_index(int)
        :return: List of parents transferred
        """
        pass

    def overlap(self, line_data_a, line_data_b):
        return line_data_a['seqid'] == line_data_b['seqid'] and (line_data_a['start'] <= line_data_b['start'] and line_data_b['start'] <= line_data_a['end'] or
                line_data_a['start'] <= line_data_b['end'] and line_data_b['end'] <= line_data_a['end'] or
                line_data_b['start'] <= line_data_a['start'] and line_data_a['end'] <= line_data_b['end'])

    def remove(self, line_data, root_type=None):
        """
        Marks line_data and all of its associated feature's 'line_status' as 'removed', does not actually remove the line_data from the data structure.
        The write function checks the 'line_status' when writing the gff file.
        Find the root parent of line_data of type root_type, remove all of its descendants.
        If the root parent has a parent with no children after the remove, remove the root parent's parent recursively.

        :param line_data:
        :param root_type:
        :return:
        """
        roots = [ld for ld in self.ancestors(line_data) if (root_type and ld['line_type'] == root_type) or (not root_type and not ld['parents'])] or [line_data]
        for root in roots:
            root['line_status'] = 'removed'
            root_descendants = self.descendants(root)
            for root_descendant in root_descendants:
                root_descendant['line_status'] = 'removed'
            root_ancestors = self.ancestors(root) # BFS, so we will process closer ancestors first
            for root_ancestor in root_ancestors:
                if len([ld for ld in root_ancestor['children'] if ld['line_status'] != 'removed']) == 0: # if all children of a root_ancestor is removed
                    # remove this root_ancestor
                    root_ancestor['line_status'] = 'removed'


    def fix(self):
        pass

    def write(self, gff_file, embed_fasta=None, fasta_char_limit=None):
        gff_fp = gff_file
        if isinstance(gff_file, str):
            gff_fp = open(gff_file, 'w')

        wrote_sequence_region = set()
        # build sequence region data
        sequence_regions = {}
        if self.fasta_external:
            for seqid in self.fasta_external:
                sequence_regions[seqid] = (1, len(self.fasta_external[seqid]['seq']))
        elif self.fasta_embedded:
            for seqid in self.fasta_embedded:
                sequence_regions[seqid] = (1, len(self.fasta_embedded[seqid]['seq']))
        else:
            pass

        wrote_lines = set()
        field_keys = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']
        reserved_attributes = ['ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from', 'Note', 'Dbxref', 'Ontology_term', 'Is_circular']
        attributes_sort_map = defaultdict(int, zip(reserved_attributes, range(len(reserved_attributes), 0, -1)))
        def write_feature(line_data):
            if line_data['line_status'] == 'removed':
                return
            field_list = [str(line_data[k]) for k in field_keys]
            attribute_list = []
            try:
                for k, v in sorted(line_data['attributes'].items(), key=lambda x: attributes_sort_map[x[0]], reverse=True):
                    if isinstance(v, list):
                        v = ','.join(v)
                    attribute_list.append('%s=%s' % (str(k), str(v)))
                field_list.append(';'.join(attribute_list))
                gff_fp.write('\t'.join(field_list) + '\n')
                wrote_lines.add(line_data['line_index'])
            except:
                logger.warning('[Missing Attributes] Program failed.\n\t\t- Line {0:s}: {1:s}'.format(str(line_data['line_index']+1), line_data['line_raw']))
        # write directives
        ignore_directives = ['##sequence-region', '###', '##FASTA']
        directives_lines = [line_data for line_data in self.lines if line_data['line_type'] == 'directive' and line_data['directive'] not in ignore_directives]
        for directives_line in directives_lines:
            gff_fp.write(directives_line['line_raw'])

        # write features
        # get a list of root nodes
        root_lines = [line_data for line_data in self.lines if line_data['line_type'] == 'feature' and not line_data['parents']]

        for root_line in root_lines:
            lines_wrote = len(wrote_lines)
            if root_line['line_index'] in wrote_lines:
                continue
            # write #sequence-region if new seqid
            if root_line['seqid'] not in wrote_sequence_region:
                if root_line['seqid'] in sequence_regions:
                    gff_fp.write('##sequence-region %s %d %d\n' % (root_line['seqid'], sequence_regions[root_line['seqid']][0], sequence_regions[root_line['seqid']][1]))
                wrote_sequence_region.add(root_line['seqid'])
            try:
                root_feature = self.features[root_line['attributes']['ID']]
            except KeyError:
                root_feature = [root_line]
            for line_data in root_feature:
                write_feature(line_data)
            descendants = self.descendants(root_line)
            for descendant in descendants:
                if descendant['line_index'] in wrote_lines:
                    continue
                write_feature(descendant)
            # check if we actually wrote something
            if lines_wrote != len(wrote_lines):
                gff_fp.write('###\n')
        # write fasta
        fasta = embed_fasta or self.fasta_external or self.fasta_embedded
        if fasta and embed_fasta != False:
            gff_fp.write('##FASTA\n')
            fasta_dict_to_file(fasta, gff_fp, line_char_limit=fasta_char_limit)

        if isinstance(gff_file, str):
            gff_fp.close()

    def sequence(self, line_data, child_type=None, reference=None):
        """
        Get the sequence of line_data, according to the columns 'seqid', 'start', 'end', 'strand'.
        Requires fasta reference.
        When used on 'mRNA' type line_data, child_type can be used to specify which kind of sequence to return:
        * child_type=None:  pre-mRNA, returns the sequence of line_data from start to end, reverse complement according to strand. (default)
        * child_type='exon':  mature mRNA, concatenates the sequences of children type 'exon'.
        * child_type='CDS':  coding sequence, concatenates the sequences of children type 'CDS'. Use the helper
                             function translate(seq) on the returned value to obtain the protein sequence.

        :param line_data: line_data(dict) with line_data['line_index'] or line_index(int)
        :param child_type: None or feature type(string)
        :param reference: If None, will use self.fasta_external or self.fasta_embedded(dict)
        :return: sequence(string)
        """
        # get start node
        reference = reference or self.fasta_external or self.fasta_embedded
        if not reference:
            raise Exception('External or embedded fasta reference needed')
        try:
            line_index = line_data['line_index']
        except TypeError:
            line_index = self.lines[line_data]['line_index']
        ld = self.lines[line_index]
        if ld['line_type'] != 'feature':
            return None
        seq = reference[ld['seqid']][ld['start']-1:ld['end']]
        if ld['strand'] == '-':
            seq = complement(seq[::-1])
        return seq

    def type_tree(self):
        class node(object):
            def __init__(self, value, children=None):
                self.value = value or ''
                self.children = children or set()

            def __repr__(self, level=0):
                ret = '\t' * level + repr(self.value) + '\n'
                for child in sorted(list(self.children), key=lambda x: x.value):
                    ret += child.__repr__(level+1)
                return ret
        root_set = set()
        node_dict = {}
        feature_line_list = [line_data for line_data in self.lines if line_data['line_type'] == 'feature']
        for line_data in feature_line_list:
            if len(line_data['children']) > 0:
                parent_type = line_data['type']
                if parent_type not in node_dict:
                    node_dict[parent_type] = node(parent_type)
                if len(line_data['parents']) == 0:
                    root_set.add(node_dict[parent_type])
                for child_ld in line_data['children']:
                    child_type = child_ld['type']
                    if child_type not in node_dict:
                        node_dict[child_type] = node(child_type)
                    if parent_type == child_type and child_type == 'mRNA':
                        print(line_data['line_index'], child_ld['line_index'])
                    else:
                        node_dict[parent_type].children.add(node_dict[child_type])
        return sorted(list(root_set), key=lambda x: x.value)

try:
    from collections import OrderedDict
except ImportError:
    # Backport of OrderedDict() class that runs on Python 2.4, 2.5, 2.6, 2.7 and pypy.
    # Passes Python2.7's test suite and incorporates all the latest updates.

    try:
        from thread import get_ident as _get_ident
    except ImportError:
        from dummy_thread import get_ident as _get_ident

    try:
        from _abcoll import KeysView, ValuesView, ItemsView
    except ImportError:
        pass


    class OrderedDict(dict):
        'Dictionary that remembers insertion order'
        # An inherited dict maps keys to values.
        # The inherited dict provides __getitem__, __len__, __contains__, and get.
        # The remaining methods are order-aware.
        # Big-O running times for all methods are the same as for regular dictionaries.

        # The internal self.__map dictionary maps keys to links in a doubly linked list.
        # The circular doubly linked list starts and ends with a sentinel element.
        # The sentinel element never gets deleted (this simplifies the algorithm).
        # Each link is stored as a list of length three:  [PREV, NEXT, KEY].

        def __init__(self, *args, **kwds):
            '''Initialize an ordered dictionary.  Signature is the same as for
            regular dictionaries, but keyword arguments are not recommended
            because their insertion order is arbitrary.

            '''
            if len(args) > 1:
                raise TypeError('expected at most 1 arguments, got %d' % len(args))
            try:
                self.__root
            except AttributeError:
                self.__root = root = []                     # sentinel node
                root[:] = [root, root, None]
                self.__map = {}
            self.__update(*args, **kwds)

        def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
            'od.__setitem__(i, y) <==> od[i]=y'
            # Setting a new item creates a new link which goes at the end of the linked
            # list, and the inherited dictionary is updated with the new key/value pair.
            if key not in self:
                root = self.__root
                last = root[0]
                last[1] = root[0] = self.__map[key] = [last, root, key]
            dict_setitem(self, key, value)

        def __delitem__(self, key, dict_delitem=dict.__delitem__):
            'od.__delitem__(y) <==> del od[y]'
            # Deleting an existing item uses self.__map to find the link which is
            # then removed by updating the links in the predecessor and successor nodes.
            dict_delitem(self, key)
            link_prev, link_next, key = self.__map.pop(key)
            link_prev[1] = link_next
            link_next[0] = link_prev

        def __iter__(self):
            'od.__iter__() <==> iter(od)'
            root = self.__root
            curr = root[1]
            while curr is not root:
                yield curr[2]
                curr = curr[1]

        def __reversed__(self):
            'od.__reversed__() <==> reversed(od)'
            root = self.__root
            curr = root[0]
            while curr is not root:
                yield curr[2]
                curr = curr[0]

        def clear(self):
            'od.clear() -> None.  Remove all items from od.'
            try:
                for node in self.__map.itervalues():
                    del node[:]
                root = self.__root
                root[:] = [root, root, None]
                self.__map.clear()
            except AttributeError:
                pass
            dict.clear(self)

        def popitem(self, last=True):
            '''od.popitem() -> (k, v), return and remove a (key, value) pair.
            Pairs are returned in LIFO order if last is true or FIFO order if false.

            '''
            if not self:
                raise KeyError('dictionary is empty')
            root = self.__root
            if last:
                link = root[0]
                link_prev = link[0]
                link_prev[1] = root
                root[0] = link_prev
            else:
                link = root[1]
                link_next = link[1]
                root[1] = link_next
                link_next[0] = root
            key = link[2]
            del self.__map[key]
            value = dict.pop(self, key)
            return key, value

        # -- the following methods do not depend on the internal structure --

        def keys(self):
            'od.keys() -> list of keys in od'
            return list(self)

        def values(self):
            'od.values() -> list of values in od'
            return [self[key] for key in self]

        def items(self):
            'od.items() -> list of (key, value) pairs in od'
            return [(key, self[key]) for key in self]

        def iterkeys(self):
            'od.iterkeys() -> an iterator over the keys in od'
            return iter(self)

        def itervalues(self):
            'od.itervalues -> an iterator over the values in od'
            for k in self:
                yield self[k]

        def iteritems(self):
            'od.iteritems -> an iterator over the (key, value) items in od'
            for k in self:
                yield (k, self[k])

        def update(*args, **kwds):
            '''od.update(E, **F) -> None.  Update od from dict/iterable E and F.

            If E is a dict instance, does:           for k in E: od[k] = E[k]
            If E has a .keys() method, does:         for k in E.keys(): od[k] = E[k]
            Or if E is an iterable of items, does:   for k, v in E: od[k] = v
            In either case, this is followed by:     for k, v in F.items(): od[k] = v

            '''
            if len(args) > 2:
                raise TypeError('update() takes at most 2 positional '
                                'arguments (%d given)' % (len(args),))
            elif not args:
                raise TypeError('update() takes at least 1 argument (0 given)')
            self = args[0]
            # Make progressively weaker assumptions about "other"
            other = ()
            if len(args) == 2:
                other = args[1]
            if isinstance(other, dict):
                for key in other:
                    self[key] = other[key]
            elif hasattr(other, 'keys'):
                for key in other.keys():
                    self[key] = other[key]
            else:
                for key, value in other:
                    self[key] = value
            for key, value in kwds.items():
                self[key] = value

        __update = update  # let subclasses override update without breaking __init__

        __marker = object()

        def pop(self, key, default=__marker):
            '''od.pop(k[,d]) -> v, remove specified key and return the corresponding value.
            If key is not found, d is returned if given, otherwise KeyError is raised.

            '''
            if key in self:
                result = self[key]
                del self[key]
                return result
            if default is self.__marker:
                raise KeyError(key)
            return default

        def setdefault(self, key, default=None):
            'od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od'
            if key in self:
                return self[key]
            self[key] = default
            return default

        def __repr__(self, _repr_running={}):
            'od.__repr__() <==> repr(od)'
            call_key = id(self), _get_ident()
            if call_key in _repr_running:
                return '...'
            _repr_running[call_key] = 1
            try:
                if not self:
                    return '%s()' % (self.__class__.__name__,)
                return '%s(%r)' % (self.__class__.__name__, self.items())
            finally:
                del _repr_running[call_key]

        def __reduce__(self):
            'Return state information for pickling'
            items = [[k, self[k]] for k in self]
            inst_dict = vars(self).copy()
            for k in vars(OrderedDict()):
                inst_dict.pop(k, None)
            if inst_dict:
                return (self.__class__, (items,), inst_dict)
            return self.__class__, (items,)

        def copy(self):
            'od.copy() -> a shallow copy of od'
            return self.__class__(self)

        @classmethod
        def fromkeys(cls, iterable, value=None):
            '''OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S
            and values equal to v (which defaults to None).

            '''
            d = cls()
            for key in iterable:
                d[key] = value
            return d

        def __eq__(self, other):
            '''od.__eq__(y) <==> od==y.  Comparison to another OD is order-sensitive
            while comparison to a regular mapping is order-insensitive.

            '''
            if isinstance(other, OrderedDict):
                return len(self)==len(other) and self.items() == other.items()
            return dict.__eq__(self, other)

        def __ne__(self, other):
            return not self == other

        # -- the following methods are only used in Python 2.7 --

        def viewkeys(self):
            "od.viewkeys() -> a set-like object providing a view on od's keys"
            return KeysView(self)

        def viewvalues(self):
            "od.viewvalues() -> an object providing a view on od's values"
            return ValuesView(self)

        def viewitems(self):
            "od.viewitems() -> a set-like object providing a view on od's items"
            return ItemsView(self)
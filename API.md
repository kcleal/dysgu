
# API


### _class_ dysgu.DysguSV(ref_genome, bam, sample_name='sample', \*\*kwargs)
This class is the main interface for calling structural variants using dysgu. To initialize DysguSV,
provide a reference genome and bam file, via the pysam library. It is also recommended to provide a sample_name
which will be used when saving results to a vcf file.


* **Parameters**

    
    * **ref_genome** (*pysam.FastaFile*) – A reference genome object from the pysam library


    * **bam** (*pysam.AlignmentFile*) – An alignment file from the pysam library


    * **sample_name** (*str*) – The sample name to use in the vcf file


    * **kwargs** (*dict*) – Key-word arguments to modify default options of dysgu


```python
import pysam
from dysgu import DysguSV

# open a reference genome and alignment file using pysam
bam = pysam.AlignmentFile('test.bam', 'rb')
genome = pysam.FastaFile('ref.fasta')

# Initialise dysgu
dysgu = DysguSV(genome, bam)

# Call SVs at a genomic location
df = dysgu("chr1:1-1000000")
```

To options can be provided during initialisation as key-word arguments:

```python
dysgu = DysguSV(genome, bam, min_support=5, mq=20)
```


#### \__call__(region, sort_df=True)
Call SVs using dysgu


* **Parameters**

    
    * **region** (*str** or **pysam Iterator*) – The genomic region to call SVs from


    * **sort_df** (*bool*) – Sort the retured dataframe



* **Returns**

    Dataframe of called SVs



* **Return type**

    pandas.DataFrame


```python
dysgu = DysguSV(ref, bam)
df = dysgu("chr1:10000-50000")

# Using a pysam iterator
df = dysgu(bam.fetch("chr1", 0, 500000))
```


#### apply_model(df)
Apply a machine leaning model to the dataframe. The model configuration is determined by the options set on
the DysguSV class. For example, to use a non diploid model, first set diploid=False:


* **Parameters**

    **df** (*pandas.DataFrame*) – The input dataframe to apply the machine learning model to



* **Returns**

    Dataframe with a modified ‘prob’ column



* **Return type**

    pandas.DataFrame


```python
dysgu = DysguSV(ref, bam)
dysgu.set_option("diploid", False)
df = dysgu("chr1:10000-50000")
```


#### call_bed_regions(regions)
Call SVs from a list of bed regions. Note,
bed regions should be sorted by genome starting position, and be non-overlapping. To create a suitable input
for this function see dysgu.load_bed() and dysgu.merge_intervals() functions.


* **Parameters**

    **regions** (*iterable*) – An iterable of bed regions



* **Returns**

    Dataframe of called SVs



* **Return type**

    pandas.DataFrame or None


```python
from dysgu import load_bed, merge_intervals

# load, merge and sort intervals
bed = load_bed('test.bed')
bed = merge_intervals(bed, srt=True)

dysgu = DysguSV(ref, bam)
df = dysgu.call_bed_regions(bed)
```


#### set_option(option, value=None)
Change option(s) for dysgu.


* **Parameters**

    
    * **option** (*str*) – The name of the option


    * **value** (*object*) – The value of the option



* **Returns**

    None



* **Return type**

    None


```python
dysgu.set_option("min_support", 10)

# Or provide a mapping of arguments:
dysgu.set_option({"min_support": 10, "mq": 20, "min_size": 100})
```


#### to_vcf(dataframe, output_file)
Save dysgu SV calls to a vcf file


* **Parameters**

    
    * **dataframe** (*pandas.DataFrame*) – A dataframe of called SVs from dysgu


    * **output_file** (*file*) – The file handle to write the vcf file to



* **Returns**

    None



* **Return type**

    None


```python
with open(path, "w") as out:
    dysgu.to_vcf(passed, out)
```


### dysgu.dysgu_default_args()
Returns the default arguments used by dysgu
:return: A dict of available arguments
:rtype: dict


### dysgu.load_dysgu_vcf(path, drop_na_columns=True)
Load a vcf file from dysgu


* **Parameters**

    
    * **path** (*str*) – The path to the vcf file


    * **drop_na_columns** (*bool*) – Drop columns that are all NAN



* **Returns**

    A dataframe of SVs



* **Return type**

    pandas.DataFrame



### dysgu.merge_dysgu_df(\*dataframes, merge_distance=500, pick_best=True, add_partners=True)
Merge calls from dysgu. Input is one or more dataframes with dysgu calls.


* **Parameters**

    
    * **dataframes** (*pandas.DataFrame*) – The input dataframes of dysgu calls to merge


    * **merge_distance** (*int*) – The merging distance, SVs closer than this spacing will be candidates for merging


    * **pick_best** (*bool*) – A single best SV is chosen for each cluster


    * **add_partners** (*bool*) – Add information to the output detailing which SVs were merged



* **Returns**

    The merged data



* **Return type**

    pandas.DataFrame



### dysgu.merge_intervals()
Merge a list of intervals, the expected format is a 3-tuple e.g. (chromosome, start, end). If add_indexes is
set to True, merge_intervals expects a 4-tuple with the last item corresponding to an index variable


* **Parameters**

    
    * **intervals** (*iterable*) – The list of intervals to merge


    * **srt** (*bool*) – Sort the intervals by chromosome and start position


    * **pad** (*int*) – Add a padding to intervals before merging. E.g. pad=10 subtracts 10 from start and adds 10 to end of interval


    * **add_indexes** (*bool*) – Add the indexes of merged intervals to the output



* **Returns**

    list of merged intervals



* **Return type**

    list


```python
>>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
>>> [['chr1', 1, 5], ['chr2', 3, 5]]
```

```python
>>> a = [("chr1", 1, 10, 0), ("chr1", 9, 11, 1), ("chr1", 20, 30, 2)]
>>> merge_intervals(a, add_indexes=True)
>>> [('chr1', 1, 11, [0, 1]), ['chr1', 20, 30, [2]]]
```

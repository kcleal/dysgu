
# FROM quay.io/pypa/manylinux2014_x86_64
FROM python:3.9

MAINTAINER Kez Cleal clealk@cardiff.ac.uk

USER root

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --upgrade pip; pip install numpy; pip install dysgu --upgrade

CMD ["/bin/sh"]


# Docker user guide
# -----------------
## Make an 'input_folder', put your bam file and reference genome inside:

# ./input_folder/
#     sample.bam
#     reference.fasta
#     reference.fasta.fai

## Make an 'output_folder'
# mkdir output_folder

## Set up docker
# docker pull kcleal/dysgu
# docker run -it \
#      --mount src="${PWD}"/input_folder,target=/input_folder,type=bind \
#      --mount src="${PWD}"/output_folder,target=/output_folder,type=bind \
#      kcleal/dysgu
# cd input_folder

## Run dysgu:

# dysgu run reference.fasta wd sample.bam > ../output_folder/sample.vcf
# exit
# sudo docker build -t kcleal/dysgu .
FROM --platform=linux/amd64 python:3.11-bullseye

MAINTAINER Kez Cleal clealk@cardiff.ac.uk

USER root

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install dysgu and required Python dependencies
RUN pip install --upgrade pip && \
    pip install "superintervals>=0.2.10" && \
    pip install dysgu

# Install samtools 1.18 for connveience
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xvf samtools-1.18.tar.bz2 && \
    rm samtools-1.18.tar.bz2 && \
    mv samtools-1.18 samtools && \
    cd samtools && \
    ./configure && \
    make -j$(nproc) && \
    make install && \
    cd ../

# Install bcftools 1.18 for connveience
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && \
    tar -xvf bcftools-1.18.tar.bz2 && \
    rm bcftools-1.18.tar.bz2 && \
    mv bcftools-1.18 bcftools && \
    cd bcftools && \
    ./configure && \
    make -j$(nproc) && \
    make install && \
    cd ../

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

#!/bin/bash

OPTIND=1
config_args=""
threads="1"
htslib_folder="./dysgu/htslib"
while getopts "h?j:l:" opt; do
    case "$opt" in
    h|\?)
        echo "Options: -j number of build threads; -l path to external htslib folder"
        exit 0
        ;;
    j)  threads=$OPTARG
        ;;
    l)  htslib_folder=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

echo "Build threads:" $threads
echo "htslib_folder:" $htslib_folder
echo "Extra configure args:" $@

cd dysgu
wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
tar -xvf htslib-1.18.tar.bz2 && rm htslib-1.18.tar.bz2 && mv htslib-1.18 htslib
cd ../

if [[ $htslib_folder == "./dysgu/htslib" ]]
then
  echo "Building htslib"
  cd ./dysgu/htslib
  autoreconf -i
  ./configure
#  autoheader
#  autoconf
#  ./configure $@
  make -j$threads
  cd ../../
fi

echo "Installing dependencies"
pip3 --version
pip3 install -r requirements.txt

echo "Installing dysgu"
python3 --version
python3 setup.py install

dysgu --version

echo "Done"

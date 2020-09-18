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


if [[ $htslib_folder == "./dysgu/htslib" ]]
then
  echo "Building htslib"
  cd ./dysgu/htslib
  autoheader
  autoconf
  ./configure $@
  make -j$threads
  cd ../../
fi

echo "Installing dependencies"
pip install -r requirements.txt

echo "Installing dysgu"
python setup.py install

dysgu --version
#dysgu test
echo "Done"

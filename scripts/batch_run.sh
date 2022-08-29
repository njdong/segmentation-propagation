#!/bin/bash


prefix="[Batch Propagation]"
version="1.0.0"
section_marker="$prefix =================================================="

echo $section_marker
echo "$prefix   Version: $version"
echo $section_marker

# validation
if test -f $1; then
  echo "$prefix FileList: $1"
else
  echo "$prefix FileList does not exist: $1"
  exit -1
fi

if test -d $2; then 
  echo "$prefix Propagation project DIR: $2"
else
  echo "$prefix Propagation project DIR does not exist: $2"
  exit -1
fi

if test -f $1; then
  echo "$prefix Configuration file: $3"
else
  echo "$prefix Configuration file does not exist: $3"
  exit -1
fi

if test -d $4; then 
  echo "$prefix Output Location: $4"
else
  echo "$prefix Output Dir does not exist. Creating one..."
  mkdir $4
fi


# use more meaningful variable names for the arguments
filelist=$1
propagation_dir=$2
config=$3
outdir=$4

# parse the file list and run foreach tag
while IFS="," read -r tag img seg tpr tpt
do
  echo $section_marker
  echo "$prefix    Starting run for tag: $tag"
  echo $section_marker
  echo "$prefix -- Image: $img"
  echo "$prefix -- Segmentation: $seg"
  echo "$prefix -- Reference TP: $tpr"
  echo "$prefix -- Target TPs: $tpt"

  # create a tag specific output folder
  tagout="$outdir/$tag"
  if test -d $tagout; then 
    :
  else
    mkdir $tagout
  fi

  # run the python script
  log="$tagout/$tag.log"
  echo "$prefix -- Start running propagation..."
  echo "$prefix -- Output is redirected to $log"

  PYTHONPATH="$propagation_dir/src" \
  python3 "$propagation_dir/propagation.py" \
  $img \
  $seg \
  $tag \
  $tpr \
  $tpt \
  "$tagout" \
  "$config" &> $log

  echo "$prefix -- Propagation completed! Return code: $?"
  
done < $filelist
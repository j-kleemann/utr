#!/bin/bash
if [ -z $TMUX ]; then
  echo "Script is not run within a tmux session! Exiting..."
  exit
fi
echo Preparing to loop-simulate an angular distribution efficiency
cd "$(dirname "$0")/.."
cmake -DUSE_GPS=OFF -DDETECTOR_CONSTRUCTION=150Nd_152Sm_775_784     .  || exit
make -j$(nproc --all) || exit

outputdir="outputdir"
threads=$(nproc --all)
for sim in "Macros/150Nd_Ang_2142_1+_2+_Loop.mac" "Macros/150Nd_Ang_2864_1-_2+_Loop.mac" "Macros/150Nd_Ang_2864_1+_2+_Loop.mac" ; do
  for file in "$outputdir"/utr*_t*.root; do
    if [ -e "$file" ] ; then 
      echo "a $file file already exists in output directory \"$outputdir\"! Exiting..."
      exit
    fi
  done

  echo "Starting simulation of $sim"
  if [[ $(uname -a) == *"Microsoft"*  ]]; then 
    ./utr -t $threads -m $sim -o $outputdir
  else 
    nice -n 18 ./utr -t $threads -m $sim -o $outputdir
  fi
  echo "Simulation $sim ended"

  for file in "$outputdir"/utr*_t*.root; do
    if [ -e "$outputdir/$(basename $sim | sed -e "s|.mac||g")_run$(basename $file| sed -e "s|utr||g")" ]; then
      echo "\"$outputdir/$(basename $sim | sed -e "s|.mac||g")_run$(basename $file| sed -e "s|utr||g")\" already exits! Skipping renaming of \"$(basename $file)\""
    else
      mv "$file" "$outputdir/$(basename $sim | sed -e "s|.mac||g")_run$(basename $file| sed -e "s|utr||g")"
    fi
  done
done

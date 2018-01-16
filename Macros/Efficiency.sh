#!/bin/bash
if [ -z $TMUX ]; then
  echo "Script is not run within a tmux session! Exiting..."
  exit
fi
echo Preparing to simulate an isotropic efficiency
cd "$(dirname "$0")/.."
cmake -DUSE_GPS=ON  -DDETECTOR_CONSTRUCTION=150Nd_152Sm_775_784_Eff .  || exit
make -j$(nproc --all) || exit

outputdir="outputdir"
threads=$(nproc --all)
for sim in "Macros/150Nd_Eff_2994_keV.mac" "Macros/150Nd_Eff_2864_keV.mac" "Macros/150Nd_Eff_2318_keV.mac" "Macros/150Nd_Eff_2142_keV.mac" ; do
  for i in $(seq 0 $(($threads-1))) ; do
    if [ -e "$outputdir/utr0_t$i.root" ]; then
      echo "a utr0_t$i.root file already exists in output directory \"$outputdir\"! Exiting..."
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

  for i in $(seq 0 $(($threads-1))) ; do
    if [ -e "$outputdir/$(basename $sim | sed -e "s|.mac||g")_t$i.root" ]; then
      echo "\"$outputdir/$(basename $sim | sed -e "s|.mac||g")_t$i.root\" already exits! Skipping renaming of \"utr0_t$i.root\""
    else
      mv "$outputdir/utr0_t$i.root" "$outputdir/$(basename $sim | sed -e "s|.mac||g")_t$i.root"
    fi
  done
done

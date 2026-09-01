#!/bin/bash


for (( i =1; i<=3; i++)); do

  for ((j=1; j<=3; j++)); do

    for ((k=0;k<=3;k++)); do

#      Rscript getQuantileScale.R $i $j $k
      sbatch  --wrap="R --no-save --no-restore < getQuantileScale.R --args $i $j $k" --output=outScale.txt

    done

  done

done


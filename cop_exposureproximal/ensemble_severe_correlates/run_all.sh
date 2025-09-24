#!/bin/bash



for (( i =1; i<=3; i++)); do

  for ((j=1; j<=3; j++)); do

    for ((k=0;k<=3;k++)); do

      #Rscript estimated3\_Scale\_LRT2.R $i $j $k

      sbatch  --wrap="R --no-save --no-restore < estimated3_Scale_LRT2.R --args $i $j $k" --output=outScale.txt

    done

  done

done



for (( i =1; i<=3; i++)); do

  for ((j=1; j<=3; j++)); do

    for ((k=0;k<=3;k++)); do

      for ((l=1;l<=500;l++));do

        sbatch  --wrap="R --no-save --no-restore < estimated3_Scale_boot_LRT2.R --args $i $j $k $l" --output=outScale.txt

      done

    done

  done

done



Rscript ComputeSimVE_Scale.R

Rscript PlotFig5Revision.R




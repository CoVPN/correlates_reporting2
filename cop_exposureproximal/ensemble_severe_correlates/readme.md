\# ENSEMBLE severe COVID exposure-proximal correlate analysis





\## Usage





\#!/bin/bash



for (( i =1; i<=3; i++)); do

&nbsp;   for ((j=1; j<=3; j++)); do

&nbsp;       for ((k=0;k<=3;k++)); do

&nbsp;	Rscript getinitial.R $i $j $k

&nbsp;              sbatch  --wrap="R --no-save --no-restore < getinitial.R --args $i $j $k" --output=outScale.txt

&nbsp;       done

&nbsp;   done

done



for (( i =1; i<=3; i++)); do

&nbsp;   for ((j=1; j<=3; j++)); do

&nbsp;       for ((k=0;k<=3;k++)); do

&nbsp;            Rscript estimated3\_Scale\_LRT2.R $i $j $k

&nbsp;           #sbatch  --wrap="R --no-save --no-restore < estimated3\_Scale\_LRT2.R --args $i $j $k" --output=outScale.txt

&nbsp;       done

&nbsp;   done

done



for (( i =1; i<=3; i++)); do

&nbsp;   for ((j=1; j<=3; j++)); do

&nbsp;       for ((k=0;k<=3;k++)); do

&nbsp;           for ((l=1;l<=500;l++));do

&nbsp;               sbatch  --wrap="R --no-save --no-restore < estimated3\_Scale\_boot\_LRT2.R --args $i $j $k $l" --output=outScale.txt

&nbsp;           done

&nbsp;       done

&nbsp;   done

done



sbatch  --wrap="R --no-save --no-restore < ComputeSimVE\_Scale.R" --output=outScale.txt

sbatch  --wrap="R --no-save --no-restore < PlotFig5Revision.R" --output=outScale.txt




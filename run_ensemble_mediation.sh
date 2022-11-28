for TRIAL in janssen_pooled_partA janssen_na_partA janssen_la_partA janssen_sa_partA
do
    sbatch --export=TRIAL=${TRIAL} --wrap "make -k -C cop_mediation all"
done

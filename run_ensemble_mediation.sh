for TRIAL in janssen_pooled_partA  janssen_la_partA
do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd cop_mediation && Rscript code/cop_mediation.R D29IncludeNotMolecConfirmed"
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd cop_mediation && Rscript code/cop_mediation.R D29SevereIncludeNotMolecConfirmed"
done

for TRIAL in janssen_pooled_partAsenior janssen_pooled_partAnonsenior janssen_na_partA janssen_na_partAsenior janssen_na_partAnonsenior janssen_la_partAsenior janssen_la_partAnonsenior janssen_sa_partA janssen_sa_partAnonsenior
do 
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd cop_mediation && Rscript code/cop_mediation.R D29IncludeNotMolecConfirmed"
done

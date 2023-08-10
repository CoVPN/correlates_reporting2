    Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
    Rscript code/cor_coxph.R D29SevereIncludeNotMolecConfirmed
    Rscript code/cor_coxph.R D29ModerateIncludeNotMolecConfirmed

    janssen_pooled_partA  janssen_na_partA janssen_la_partA janssen_sa_partA 

for run_survtmle in TRUE FALSE
  do
  for impute_placebo in TRUE FALSE
    do
    for outcome in D29IncludeNotMolecConfirmed D29SevereIncludeNotMolecConfirmed D29ModerateIncludeNotMolecConfirmed
      do
      # POOLED
      for TRIAL in janssen_pooled_partA
        do
        for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50
          do
          sbatch --export=TRIAL=${TRIAL} \
          --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R ${outcome} ${run_survtmle} ${impute_placebo} ${MARKER}"
        done
      done
      # LATIN AMERICA
      for TRIAL in janssen_la_partA
        do
        for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50la
          do
          sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${run_survtmle} ${impute_placebo} ${MARKER}"
        done
      done
      # SOUTH AFRICA
      for TRIAL in janssen_sa_partA
        do
        for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50sa
          do
          sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${run_survtmle} ${impute_placebo} ${MARKER}"
        done
      done
      # NORTH AMERICA 
      for TRIAL in janssen_na_partA
        do
        for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50
          do
          sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${run_survtmle} ${impute_placebo} ${MARKER}"
        done
      done
    done # END OUTCOME
  done # END IMPUTE PLACEBO
done # END SURVTMLE

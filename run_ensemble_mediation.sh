for run_survtmle in TRUE FALSE
  do
  for impute_placebo in TRUE
    do
    for truncate_tf_day in TRUE
    do
      for compatible_total_ve in FALSE
      do
      for outcome in D29IncludeNotMolecConfirmed D29SevereIncludeNotMolecConfirmed D29ModerateIncludeNotMolecConfirmed
        do
        # POOLED
        for TRIAL in janssen_pooled_partA
          do
          for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50
            do
            sbatch --export=TRIAL=${TRIAL} \
            --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R ${outcome} ${run_survtmle} ${impute_placebo} ${MARKER} ${truncate_tf_day} ${compatible_total_ve}"
          done
        done
        # LATIN AMERICA
        for TRIAL in janssen_la_partA
          do
          for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50la
            do
            sbatch --export=TRIAL=${TRIAL} \
             --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R ${outcome} ${run_survtmle} ${impute_placebo} ${MARKER} ${truncate_tf_day} ${compatible_total_ve}"
          done
        done
        # SOUTH AFRICA
        for TRIAL in janssen_sa_partA
          do
          for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50sa
            do
            sbatch --export=TRIAL=${TRIAL} \
             --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R ${outcome} ${run_survtmle} ${impute_placebo} ${MARKER} ${truncate_tf_day} ${compatible_total_ve}"
          done
        done
        # NORTH AMERICA 
        for TRIAL in janssen_na_partA
          do
          for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50
            do
            sbatch --export=TRIAL=${TRIAL} \
             --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R ${outcome} ${run_survtmle} ${impute_placebo} ${MARKER} ${truncate_tf_day} ${compatible_total_ve}"
          done
        done
      done # END OUTCOME
    done # END COMPATIBLE VE
  done # END TRUNCATE TF DAY
  done # END IMPUTE PLACEBO
done # END SURVTMLE

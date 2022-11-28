# Geog regions pooled COVID-19 endpoint: 4 markers, for Overall, Senior, and non-Senior
for TRIAL in janssen_pooled_partA janssen_pooled_partAsenior janssen_pooled_partAnonsenior
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${MARKER}"
    done
done

# Latin America COVID-19 endpoint: 5 markers, for Overall, Senior, and non-Senior
for TRIAL in janssen_la_partA janssen_la_partAsenior janssen_la_partAnonsenior
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50la Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${MARKER}"
    done
done

# South Africa COVID-19 endpoint: 5 markers, for Overall, Senior, and non-Senior
for TRIAL in janssen_sa_partA janssen_sa_partAnonsenior
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50sa Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${MARKER}"
    done
done

# U.S. COVID-19 endpoint: 4 markers, for Overall, Senior, and non-Senior
for TRIAL in janssen_na_partA janssen_na_partAsenior janssen_na_partAnonsenior  
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29IncludeNotMolecConfirmed ${MARKER}"
    done
done

# Geog regions pooled severe COVID-19 endpoint: 4 markers for Overall
for TRIAL in janssen_pooled_partA
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29SevereIncludeNotMolecConfirmed ${MARKER}"
    done
done

# Latin America severe COVID-19 endpoint: 5 markers for Overall
for TRIAL in janssen_la_partA
do
    for MARKER in Day29bindSpike Day29bindRBD Day29pseudoneutid50 Day29pseudoneutid50la Day29ADCP
    do
    sbatch --export=TRIAL=${TRIAL} \
           --wrap="cd ~/correlates_reporting2/cop_mediation && /app/software/R/4.0.4-foss-2020b/bin/Rscript code/mediation.R D29SevereIncludeNotMolecConfirmed ${MARKER}"
    done
done




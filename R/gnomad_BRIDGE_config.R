studies = list()
studies$gnomad = list("study_name" = "gnomad", "vars" = "~/scratch/gnomAD/gnomad.genomes.r2.0.1.sites.%s.txt", 'null_model' = "~/software/selection2/models/obs_exp_lm.gnomAD_NFE.30x_cov.RData", "pop_size" = 7509, "AC_column_name" = 'AC_NFE')
studies$BRIDGE = list("study_name" = "BRIDGE", "vars" = "~/scratch/BRIDGE/chr%s.tsv", 'null_model' = "~/software/selection2/models/obs_exp_lm.BRIDGE.30x_cov.RData", "pop_size" = 13049, "AC_column_name" = 'AC')


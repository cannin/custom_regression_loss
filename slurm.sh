#!/bin/bash

#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mail-user=augustin_luna@hms.harvard.edu
#SBATCH --mail-type=ALL
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

#singularity exec /n/app/singularity/containers/al323/tidy-4.0.4.sif bash -c "R -e '.libPaths(c(\"~/singularity_tmp\", .libPaths())); setwd(\"~/singularity_tmp\"); r_packages <- c(\"BiocManager\", \"numDeriv\", \"gridExtra\"); lapply(r_packages, function(x) { install.packages(x)}); BiocManager::install(\"fgsea\", ask=FALSE, update=FALSE); library(numDeriv); library(fgsea); data(examplePathways); data(exampleRanks); fr <- fgseaSimple(examplePathways, exampleRanks, nperm=10, maxSize=100); gr <- grad(sin, pi); cat(runif(1), file=\"out.txt\")'"
# 
# srun --pty -p interactive -t 0-12:00 --mem 8192 /bin/bash 
# singularity exec /n/app/singularity/containers/al323/tidy-4.0.4.sif bash

echo ${SLURM_ARRAY_TASK_ID}

singularity exec /n/app/singularity/containers/al323/tidy-4.0.4.sif bash -c "R -e '.libPaths(c(\"~/singularity_rlibs\", .libPaths())); cran_packages <- c(\"BiocManager\", \"numDeriv\", \"gridExtra\"); lapply(cran_packages, function(x) { if (!(x %in% .packages(TRUE))) install.packages(x) }); bioc_packages <- c(\"fgsea\"); lapply(bioc_packages, function(x) { if (!(x %in% .packages(TRUE))) BiocManager::install(x, ask=FALSE, update=FALSE) }); setwd(\"~/fgsea_model\"); source(\"model_modified.R\")'"
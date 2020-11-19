######### script.R ####################################################
#' R script start the RnBeads analysis of the TCGA SCKM dataset. Notably,
#' we used an SGE HPC environment, but RnBeads can also be executed
#' on a standalone workstation or on a SLURM HPC
library(RnBeads)
theme_set(theme_bw())
xml.file <- "rnbOptions.xml"
arch <- new("ClusterArchitectureSGE")
arch <- setExecutable(arch,"R","/usr/bin/R")
arch <- setExecutable(arch,"Rscript","/usr/bin/Rscript")
rnb.cr <- new("RnBClusterRun",arch)
rnb.cr <- setModuleResourceRequirements(rnb.cr,c(mem_free="120G",h_vmem="120G"),"all")
rnb.cr <- setModuleNumCores(rnb.cr,4L,"all")
rnb.cr <- setModuleNumCores(rnb.cr,2L,"exploratory")
run(rnb.cr, "rnB_TCGA_SCKM", xml.file)

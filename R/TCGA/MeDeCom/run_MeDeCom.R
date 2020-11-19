######### run_MeDeCom.R ####################################################
#' Main R script to perform the MeDeCom analysis on the TCGA SKCM data.
#' The scripts uses a recently
#' published protocol (https://www.nature.com/articles/s41596-020-0369-6)
#' to perform reference-free deconvolution. Notably,
#' we used an SGE HPC environment, but MeDeCom can also be executed
#' on a standalone workstation or on a SLURM HPC
#' Input is the output of the RnBeads run
library(DecompPipeline)

options(fftempdir="/tmp/")

rnb.set <- load.rnb.set("TCGA_SKCM_report/cluster_run/import_RnBSet/")

md.res <- start.decomp.pipeline(rnb.set,
	Ks=2:15,
	lambda.grid=c(0.01,0.001,0.0001),
	work.dir="TCGA/MeDeCom/",
	analysis.name="DecompICA_age_reduced",
	id.column="Sample_ID",
	normalization="dasen",
	filter.beads=T,
	filter.intensity=T,
	min.int.quant=0.001,
	max.int.quant=0.999,
	filter.na=TRUE,
	filter.context=T,
	filter.cross.reactive=T,
	execute.lump=T,
	remove.ICA=T,
	conf.fact.ICA=c("age_at_diagnosis","gender"),
	n.markers=5000,
	filter.coverage=T,
	filter.snp=T,
	filter.sex.chromosomes=T,
	marker.selection=c("var","hybrid","random"),
	cores=10,
	cluster.submit=T,
	cluster.Rdir="/usr/bin",
	factorviz.outputs=T,
	ica.setting=c(alpha.fact=1e-5,save.report=TRUE)
)


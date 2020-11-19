library(RnBeads)
rnb.set <- load.rnb.set("melanoma_rnb_report/cluster_run/preprocessing_RnBSet/")
ph <- pheno(rnb.set)
meth.data <- meth(rnb.set)
median.samples <- list()
for(pid in unique(ph$Patient_ID_850K)){
	print(pid)
	sel.samples <- ph$Patient_ID_850K %in% pid
	if(sum(sel.samples)==1){
		median.samples[[pid]] <- c(sample=as.character(ph$Sample_ID[sel.samples]),dist=0)
		next
	}
	meth.sel <- meth.data[,sel.samples]
	med.cpgs <- apply(meth.sel,1,median,na.rm=T)
	dist.median <- as.matrix(dist(t(data.frame(meth.sel,med.cpgs))))
	min.dist <- min(dist.median[nrow(dist.median),-ncol(dist.median)])
	med.sample <- which.min(dist.median[nrow(dist.median),-ncol(dist.median)])
	med.sample <- ph$Sample_ID[sel.samples][med.sample]
	median.samples[[pid]] <- c(sample=as.character(med.sample),dist=min.dist)
}
median.samples <- unlist(median.samples)
median.samples <- data.frame(Sample_ID=median.samples[seq(1,83,2)],Distance=median.samples[seq(2,84,2)])
write.csv(median.samples,"representative_samples.csv")
s.anno <- read.csv("sample_annotation.csv")
s.anno <- s.anno[s.anno$Sample_ID%in%median.samples$Sample_ID,]
write.csv(s.anno,"sample_annotation_reduced.csv")


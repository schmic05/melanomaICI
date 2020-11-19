library(RnBeads)
my.colors <- list(ICI_response=c("#D6604D","#2ca25f","#2ca25f"),
	LUMP_estimate=rev(c("#D6604D", "#F4A582","#FDDBC7")),
	BRAF=c(mt="indianred3",wt="cyan4"),
	NRAS=c(mt="indianred3",wt="cyan4"),
	Brain_Mets=c(no="#D6604D",yes="#2ca25f"),
	Subtype=c(ALM="#E7298A",NM="#1B9E77",CUP="#E6AB02",SSM="#D95F02",AM="#D6604D",UNK="gray80"))
my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"))

rnb.set <- load.rnb.set("melanoma_rnb_report/cluster_run/preprocessing_RnBSet/")
rnb.set <- addPheno(rnb.set,pheno(rnb.set)$'LUMP estimate',"LUMP_estimate")

for(trait in names(my.colors)){
	rnb.options(colors.category=my.colors[[trait]])
	if(trait%in%c("LUMP_estimate","ICI_response")){
		rnb.options(colors.3.gradient=my.colors[[trait]])
	}	
plot <- rnb.plot.dreduction(rnb.set,point.colors=trait)+my_theme+scale_fill_continuous(na.value="gray80")
	ggsave(paste0("PCA_",trait,".pdf"),plot)
}

library(pheatmap)
my.colors <- list(ICI_response=c(no="#D6604D",yes="#2ca25f"),
	LUMP_estimate=rev(c("#D6604D", "#F4A582","#FDDBC7")),
	BRAF=c(mt="indianred3",wt="cyan4"),
	NRAS=c(mt="indianred3",wt="cyan4"),
	Brain_Mets=c(no="#D6604D",yes="#2ca25f"),
	Subtype=c(ALM="#E7298A",NM="#1B9E77",CUP="#E6AB02",SSM="#D95F02",AM="#D6604D",UNK="gray80"))
meth.data <- meth(rnb.set)
sds <- apply(meth.data,1,sd,na.rm=T)
sel.meth <- meth.data[(order(sds,decreasing=T)[sds<0.24])[1:1000],]
ph <- pheno(rnb.set)
anno.col <- ph[,"ICI_response",drop=F]
row.names(anno.col) <- ph$Patient_ID_850K
colnames(sel.meth) <- ph$Patient_ID_850K
png("heatmap_wo_SNPs.png")
pheatmap(sel.meth,annotation_col=anno.col,annotation_color=my.colors,
	col=rev(colorRampPalette(c("#67001F", 
				"#B2182B", "#D6604D", "#F4A582",
                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                "#4393C3", "#2166AC", "#053061"))(200)))
dev.off()


my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=25,color="black"))

ph <- pheno(rnb.set)
to.plot <- ph[,c("ICI_response","LUMP_estimate")]
to.plot$ICI_response <- ifelse(to.plot$ICI_response==0,"no","yes")
plot <- ggplot(to.plot,aes(x=ICI_response,y=LUMP_estimate,fill=ICI_response))+geom_boxplot()+my_theme+
	scale_fill_manual(values=my.colors[["ICI_response"]])

diff.meth <- read.csv("melanoma_rnb_report/differential_methylation_data/diffMethTable_site_cmp1.csv")
plot <- create.densityScatter(diff.meth[,c("mean.yes","mean.no")],is.special=diff.meth$diffmeth.p.val<1e-5,add.text.cor=T)+
	my_theme+theme(legend.position="none")+xlab("ICI responder")+ylab("ICI non-responder")
ggsave("ICI_response_diffMeth.png",plot,dpi=600)

diff.meth <- read.csv("melanoma_rnb_report/differential_methylation_data/diffMethTable_site_cmp2.csv")
plot <- create.densityScatter(diff.meth[,c("mean.yes","mean.no")],is.special=diff.meth$diffmeth.p.val<1e-5,add.text.cor=T)+
	my_theme+theme(legend.position="none")+xlab("Brain metastases")+ylab("no Brain metastases")
ggsave("Brain_mets_diffMeth.png",plot,dpi=600)
#sum(diff.meth$diffmeth.p.val<1e-5,na.rm=T)   
#[1] 3

diff.meth <- read.csv("melanoma_rnb_report//differential_methylation_data/diffMethTable_site_cmp3.csv")
plot <- create.densityScatter(diff.meth[,c("mean.mt","mean.wt")],is.special=diff.meth$diffmeth.p.val<1e-5,add.text.cor=T)+
	my_theme+theme(legend.position="none")+xlab("BRAF mutant")+ylab("BRAF wildtype")
ggsave("BRAF_diffMeth.png",plot,dpi=600)
#> sum(diff.meth$diffmeth.p.val<1e-5,na.rm=T) 
#[1] 1

diff.meth <- read.csv("melanoma_rnb_report/differential_methylation_data/diffMethTable_site_cmp4.csv")
plot <- create.densityScatter(diff.meth[,c("mean.mt","mean.wt")],is.special=diff.meth$diffmeth.p.val<1e-5,add.text.cor=T)+
	my_theme+theme(legend.position="none")+xlab("NRAS mutant")+ylab("NRAS wildtype")
ggsave("NRAS_diffMeth.png",plot,dpi=600)


######### plotting_scripts.R ####################################################
#' R script to plot various aspects of the MeDeCom analysis. Thus, the input to
#' the script is the output of run_MeDeCom.R
library(pheatmap)
library(MeDeCom)
library(RColorBrewer)
factorviz.output <- "DecompICA_age/FactorViz_outputs/"
K <- 6
lambda <- 0.001
cg_subset <- 1
traits <- c("LUMP_estimate","ICI_response")
s.id.col <- "Patient_ID_850K"
plot.path <- "MeDeCom/"
plot.name <- "heatmap_variable_5000_new.pdf"
max.val <- NULL
min.val <- NULL

load(paste0(factorviz.output,"medecom_set.RData"))
contris <- as.data.frame(getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subset))
load(paste0(factorviz.output,"ann_S.RData"))
colnames(contris) <- ann.S[,s.id.col]
if(is.null(max.val)){
  max.val <- max(apply(contris,2,max))
}
if(is.null(min.val)){
  min.val <- min(apply(contris,2,min))
}
breaksList <- seq(min.val,max.val,by=0.01)
sel.traits <-ann.S[,traits]
row.names(sel.traits) <- ann.S[,s.id.col]
#sel.traits$ICI_response <- ifelse(sel.traits$ICI_response==0,"no","yes")
pdf(file.path(plot.path,plot.name))
pheatmap(contris,cluster_rows = F,clustering_distance_cols = "euclidean",clustering_method="ward.D",show_colnames = T,
        breaks=breaksList,annotation_color=list(ICI_response=c(no="#D6604D",yes="#2ca25f"),LUMP_estimate=rev(c("#D6604D", "#F4A582","#FDDBC7"))),
        color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(length(breaksList)))
dev.off()


library(corrplot)
K <- 6
lambda <- 0.001
props <- getProportions(medecom.set,K=K,lambda=lambda)
sel.traits <- c("LUMP_estimate","survival")
cors <- apply(props,1,function(x){
  sapply(sel.traits,function(trait){
    trait <- ann.S[,trait]
    na.trait <- is.na(trait)
    cor(x[!na.trait],as.numeric(trait[!na.trait]))
#    cor.test(x[!na.trait],as.numeric(trait[!na.trait]))$p.value
  })
})
pdf(file.path(plot.path,"trait_association.pdf"))
corrplot(t(cors),"ellipse",addCoef.col = "black",col=rev(colorRampPalette(c("#67001F", 
				"#B2182B", "#D6604D", "#F4A582",
                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                "#4393C3", "#2166AC", "#053061"))(200)))
dev.off()

sel.traits <- c("Subtype","BRAF","Censor","ICI_response")
#ann.S$ICI_response <- factor(ann.S$ICI_response,levels=c("0","1"))
madiff <- apply(props,1,function(x){
  sapply(sel.traits,function(trait){
	print(trait)
    trait <- ann.S[,trait]
    na.trait <- is.na(trait)
    trait <- trait[!na.trait]
    x <- x[!na.trait]
#    md <- mean(x[as.character(trait)==levels(trait)[1]],na.rm=T)-mean(x[as.character(trait)==levels(trait)[2]],na.rm=T)
#	names(md) <- paste0(levels(trait)[1],"vs",levels(trait)[2])
#	md
   t.test(x[as.character(trait)==levels(trait)[1]],na.rm=T,x[as.character(trait)==levels(trait)[2]])$p.value
  })
})
to.plot <- data.frame(t(madiff),LMC=colnames(madiff))
to.plot <- melt(to.plot,id="LMC")
colnames(to.plot)[2:3] <- c("Trait","ADiffM")
mini <- min(madiff)
maxi <- max(madiff)
abs.max <- max(abs(mini),abs(maxi))
plot <- ggplot(to.plot,aes(x=Trait,y=LMC,fill=ADiffM))+geom_tile()+theme_bw()+theme(panel.grid=element_blank())+
  scale_fill_gradient2(low="#6a011f",mid="white",midpoint=0,high="#35ac35",limits=c(-abs.max,abs.max))
ggsave(file.path(plot.path,"trait_association_qual.pdf"),plot)

library(pheatmap)
library(MeDeCom)
library(RColorBrewer)
factorviz.output <- "DecompICA_age/FactorViz_outputs/"
K <- 6
lambda <- 0.001
cg_subset <- 1
traits <- c("LUMP_estimate","ICI_response")
s.id.col <- "Patient_ID_850K"
plot.path <- "MeDeCom/"
plot.name <- "heatmap_variable_5000_new.pdf"
max.val <- NULL
min.val <- NULL
my.colors <- list(ICI_response=c(no="#D6604D",yes="#2ca25f"),
	LUMP_estimate=rev(c("#D6604D", "#F4A582","#FDDBC7")),
	BRAF=c(mt="indianred3",wt="cyan4"),
	NRAS=c(mt="indianred3",wt="cyan4"),
	Brain_Mets=c(no="#D6604D",yes="#2ca25f"),
	Subtype=c(ALM="#E7298A",NM="#1B9E77",CUP="#E6AB02",SSM="#D95F02",AM="#D6604D",UNK="gray80"))
my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"))

load(paste0(factorviz.output,"medecom_set.RData"))
contris <- as.data.frame(getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subset))
load(paste0(factorviz.output,"ann_S.RData"))
colnames(contris) <- ann.S[,s.id.col]
to.plot <- data.frame(LMC4=unlist(contris[4,]),ICI=ann.S$ICI_response)
plot <- ggplot(to.plot,aes(x=ICI,y=LMC4,fill=ICI))+geom_boxplot()+my_theme+scale_fill_manual(values=my.colors[["ICI_response"]])+
	xlab("ICI response")+ylab("LMC4 proportion")


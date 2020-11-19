######### lasso_predictor.R ####################################################
#' R script to perform a pilot classification model for the LMC-based clustering
#' Input is the output directory of MeDeCom specified in run_MeDeCom.R
library(MeDeCom)
library(glmnet)
library(pROC)
# Specify the path to the MeDeCom output directory specified in run_MeDeCom.R here
load("DecompICA_age/FactorViz_outputs/meth_data.RData")
load("DecompICA_age/FactorViz_outputs/ann_S.RData")
load("DecompICA_age/FactorViz_outputs/ann_C.RData")
load("DecompICA_age/FactorViz_outputs/medecom_set.RData")

props <- getProportions(medecom.set,K=6,lambda=0.001)
output.variable <- ifelse(props[2,]>0.350360121,1,0)

lmcs <- getLMCs(medecom.set,K=6,lambda=0.001)
hypo.lmc2 <- apply(lmcs[,-2],1,median)-lmcs[,2]>0.75
hyper.lmc2 <- apply(lmcs[,-2],1,median)-lmcs[,2]<(-0.85)

input.variables <- meth.data[medecom.set@parameters$GROUP_LISTS[[1]],]
ann.C <- ann.C[medecom.set@parameters$GROUP_LISTS[[1]],]
input.variables <- input.variables[hypo.lmc2|hyper.lmc2,]
write.csv(ann.C[hyper.lmc2,],"hyper_methylated_CpGs_LMC2.csv")
write.csv(ann.C[hypo.lmc2,],"hypo_methylated_CpGs_LMC2.csv")
ann.C <- ann.C[hypo.lmc2|hyper.lmc2,]

#out.res <- list()
res.all <- c()
pred.all <- c()

for(i in 1:100){
	permu <- sample(1:ncol(input.variables),ncol(input.variables))
	used.input <- input.variables[,permu]
	used.output <- output.variable[permu]
	inner.res <- c()
	for(j in 0:9){
		if(j<9){
			sel.samples <- ((j*4)+1):((j+1)*4)
		}else{
			sel.samples <- 32:39
		}
		test.input <- t(used.input[,sel.samples])
		test.output <- used.output[sel.samples]
		train.input <- t(used.input[,-sel.samples])
		train.output <- used.output[-sel.samples]
		cv.glm <- cv.glmnet(x=train.input,y=train.output,alpha=1,family="binomial")
		mod <- glmnet(x=train.input,y=train.output,alpha=1,family="binomial",lambda=cv.glm$lambda.min)
		fit.output <- predict(mod,test.input)
		#inner.res <- c(inner.res,sum((fit.output>0)!=test.output)/length(test.output))
		res.all <- c(res.all,test.output)
		pred.all <- c(pred.all,fit.output)
	}
#	out.res[[i]] <- inner.res
}

r <- roc(res.all,pred.all)
pdf("ROC.pdf")
plot(r)
dev.off()
#mean.res <- lapply(out.res,mean)
cv.glm <- cv.glmnet(x=t(input.variables),y=output.variable,alpha=1,family="binomial")
final.model <- glmnet(x=t(input.variables),y=output.variable,alpha=1,family="binomial",lambda=cv.glm$lambda.min)

ann.C <- ann.C[which(coef(final.model)[-1]!=0),]
ann.C.gr <- makeGRangesFromDataFrame(ann.C)
genes <- unlist(rnb.get.annotation("genes","hg19"))
dists <- c()
symbs <- c()
for(i in 1:length(ann.C.gr)){
	cpg <- ann.C.gr[i]
	closest.gene <- which.min(distance(cpg,genes))
	dists <- c(dists,min(distance(cpg,genes),na.rm=T))
	symbs <- c(symbs,values(genes[closest.gene])$symbol)
}
ann.C$gene <- symbs
ann.C$distance <- dists
ann.C$coef <- coef(final.model)[-1][which(coef(final.model)[-1]!=0)]
ann.C$Chromosome <- as.character(ann.C$Chromosome)
ann.C <- rbind(c("Intercept",rep(NA,23),coef(final.model)[1]),ann.C)
write.csv(ann.C,"model_coefs.csv")


####Clustering####
setwd(" ")
#BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)
#1.2 
cellMarker <- data.table::fread("KEGG human metabolic pathway geneset.csv",data.table = F,header = T)
colnames(cellMarker)[1] <- "Metagene"
colnames(cellMarker)[2] <- "celltype"

type <- split(cellMarker,cellMarker$celltype)

cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})

save(cellMarker,file = "cellmarker_ssGSEA.Rdata")#保存中间文件

##1.3 
expr <- data.table::fread("gene_fpkm.txt",data.table = F) #读取表达文件

expr = expr %>% distinct(ID,.keep_all = T) %>% as.tibble() %>%
  column_to_rownames(var = "ID") 
expr <- as.matrix(expr)
#2. 
gsva_data <- gsva(expr,cellMarker, method = "gsva")
gsva_data<-as.data.frame(gsva_data)
write.table(gsva_data,"gsva_data.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

counts <- data.table::fread("gene_raw counts.txt",data.table = F) 
colnames(counts)[1] <- "ID"
counts= counts %>% distinct(ID,.keep_all = T) %>% as.tibble() %>%
  column_to_rownames(var = "ID") 
gene_rawcounts<-counts[,genelist]
#write.csv(gene_rawcounts,file = "gene_rawcounts_M.csv")

setwd(" ")
library(tidyverse)
library(data.table)
library(sweep)
library(dplyr)
gsva_data<-read.table("gsva_data.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gsva_data<-as.matrix(gsva_data)
rowsd=apply(gsva_data,1,sd)#sd/mad
rowsd<-as.data.frame(rowsd)
rowsd <- rowsd %>% rownames_to_column("Gene")
gsva_data<-as.data.frame(gsva_data)
gsva_data <- gsva_data %>% rownames_to_column("Gene") %>% as.data.frame()
gsva_data <- cbind(gsva_data,rowsd)

gsva_data = gsva_data %>% distinct(Gene,.keep_all = T) %>% as.tibble()

gsva_data<-arrange(gsva_data,desc(rowsd))
gsva_data<-gsva_data[-1,]
gsva_data<-gsva_data[c(1:86),-c(103:104)]#c(1:100)

rownames(gsva_data) <- NULL
gsva_data <- gsva_data %>% column_to_rownames("Gene") 
#z-score
gsva_data <- scale(gsva_data)%>% as.matrix(gsva_data)

##ConsensusClusterPlus
#install.packages("BiocManager")
#BiocManager::install('ConsensusClusterPlus')
library(tidyverse)
library(ConsensusClusterPlus)

d<-gsva_data
d1 <- sweep(d,1,apply(d,1,median))
title=("clustering") 
results = ConsensusClusterPlus(d1,
                               maxK=8,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title=title,
                               clusterAlg="km",#hc，pam，km；
                               distance="pearson",#pearson，spearman，euclidean；
                               seed=1,
                               plot="pdf")
#results[[2]][["consensusMatrix"]][1:5,1:5]
#results[[2]][["consensusTree"]]
#results[[2]][["consensusClass"]][1:5]

icl = calcICL(results,title=title,plot="pdf") 

group<-results[[3]][["consensusClass"]]
group<-as.data.frame(group)
table(group)
group$group <- factor(group$group,levels=c(1,2,3))
save(group,file = "group.Rda")


####Metabolomics####
#library(BiocManager)
#install("ropls")# 用于PLS-DA和OPLS-DA分析
library(ropls)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggthemes)
setwd("/Users/liuchuqiao/Desktop/MTC-metabolism/metabo")
data_matrix <- data.table::fread("metabolite_matrix2.csv",data.table = F,header = T)
data_matrix = data_matrix %>% distinct(V1,.keep_all = T) %>% as.tibble() %>% column_to_rownames(var = "V1")
#data_matrix <- t(data_matrix)
#write.csv(data_matrix,file = "metabolite_matrix2.csv")
data_group <- data.table::fread("metabolite_group.csv",data.table = F,header = T)
data_group <- data_group[ order (data_group$Group2,decreasing = T ), ]
identical(data_group$Sample,rownames(data_matrix))
data<-t(data_matrix[,-1])

#PCA
data_pca <- opls(t(data))
plot(data_pca,
     typeVc = "x-score",
     parAsColFcVn = groupFc,
     parPaletteVc=c("green4", "orange"))

#OPLS-DA
Group <- data_group$Group2
library(scatterplot3d)
library(dplyr)
library(tidyverse)
data_oplsda_3d <- opls(t(data), Group, predI = 1, orthoI = 2)

oplsda_scores_3d = cbind(data_oplsda_3d@scoreMN, data_oplsda_3d@orthoScoreMN)#得分矩阵
identical(rownames(oplsda_scores_3d),data_group$Sample )
oplsda_scores_3d<-as.data.frame(oplsda_scores_3d)
oplsda_scores_3d$group <- data_group$Group2 
oplsda_scores_3d = oplsda_scores_3d %>% 
  mutate(colour = case_when(oplsda_scores_3d$group == "M12" ~ '#9b9c9f',
                            oplsda_scores_3d$group == "M3" ~ '#e69a96'))
colnames(oplsda_scores_3d)[c(1,2,3)]<-c("PC1","PCo1","PCo2")


####LASSO_deg####
setwd(" ")
library(glmnet)
library(tidyverse)

dat1<-data.table::fread("gene_fpkm_ML.txt", data.table = F,header = T)
dat1 = dat1 %>% distinct(id,.keep_all = T) %>% as.tibble()%>%
  column_to_rownames(var = "id")
dat2<-data.table::fread("M3_metagene_diffsig_up.csv", data.table = F,header = T)
dat2 = dat2 %>% distinct(V1,.keep_all = T) %>% as.tibble()%>%
  column_to_rownames(var = "V1")
genelist <- intersect(rownames(dat1),rownames(dat2))
dat<-dat1[genelist,]
dat<-t(dat)
x <- as.matrix(dat)
y <- gsub("(.*)\\-(.*)\\-(.*)", "\\3", row.names(dat))

dat<-data.table::fread("surv_sig_meta_diff.csv", data.table = F,header = T)
x <- as.matrix(dat[,-c(1:3)])
y <- ifelse(dat$status == "0", 0,1)#把分组信息换成0或1

fit=glmnet(x, y, family = "binomial", alpha=1)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
plot(cvfit)
dev.off()

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
lassoGene



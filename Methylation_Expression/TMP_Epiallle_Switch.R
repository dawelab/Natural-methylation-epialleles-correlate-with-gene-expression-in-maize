##I will make use the methylation matrix generated in Epiallele Stability
library(tidyverse)
library(matrixStats)
single_pan <- data.frame(pan=df_NAM_copy[df_NAM_copy$Min==1&df_NAM_copy$Max==1,1])
single_methylation_matrix = merge(single_pan,Matrix_epiallele,by="pan",all.x = T)

#Create the expression matrix for different tissues
file_names<- list.files("/Users/x/Desktop/Data/expression")[-c(1:6,17,25,30:31,27,28,37:41,44)]
path <- "/Users/x/Desktop/Data/expression/"
for (i in 1:length(file_names))  assign(gsub(".txt","",file_names[i]), read.table(paste0(path,file_names[i]),header = T))

NAM_list <- cbind(file_names,name)
NAM_list[12,1] <- "IL14H.tpm.txt"
NAM_list[12,2] <-"Il14H"

tissue_list = names(B73.tpm)[5:14]

#Filter out the unqualified genes
Filter_expression <- function(df_class,df_tpm,tissue,NAM){ 
  df = merge(df_class,df_tpm,by="gene")[,c("pan",tissue)]
  names(df) <- c("pan",NAM)
  return(df)
}
#head(Filter_expression(B73_class,B73.tpm,"tip","B73"))

for (j in c(3:7,9,10)) {
  for (i in 1:26)  assign(paste0(NAM_list[i,2],tissue_list[j]), Filter_expression(get(paste0(NAM_list[i,2],"_class")),get(gsub(".txt","",NAM_list[i,1])),tissue_list[j],NAM_list[i,2]))  
}

embroy_NAM = NAM_list[-c(4,5,6,8,9,16),];j=1
for (i in 1:20) assign(paste0(embroy_NAM[i,2],tissue_list[j]), Filter_expression(get(paste0(embroy_NAM[i,2],"_class")),get(gsub(".txt","",embroy_NAM[i,1])),tissue_list[j],embroy_NAM[i,2]))  

tip_NAM = NAM_list[-c(22),];j=8
for (i in 1:25)  assign(paste0(tip_NAM[i,2],tissue_list[j]), Filter_expression(get(paste0(tip_NAM[i,2],"_class")),get(gsub(".txt","",tip_NAM[i,1])),tissue_list[j],tip_NAM[i,2]))  

endosperm_NAM = NAM_list[-c(4,6,9),];j=2
for (i in 1:23)  assign(paste0(endosperm_NAM[i,2],tissue_list[j]), Filter_expression(get(paste0(endosperm_NAM[i,2],"_class")),get(gsub(".txt","",endosperm_NAM[i,1])),tissue_list[j],endosperm_NAM[i,2]))  

embryo_mat <- endosperm_mat <- root_mat <- shoot_mat <- anther_mat <- base_mat <- middle_mat <- tip_mat <- ear_mat <- tassel_mat <- single_pan
for (i in 1:26) shoot_mat = merge(shoot_mat,get(paste0(NAM_list[i,2],tissue_list[4])),by="pan",all.x = T)
for (i in 1:26) anther_mat = merge(anther_mat,get(paste0(NAM_list[i,2],tissue_list[5])),by="pan",all.x = T)
for (i in 1:26) base_mat = merge(base_mat,get(paste0(NAM_list[i,2],tissue_list[6])),by="pan",all.x = T)
for (i in 1:26) middle_mat = merge(middle_mat,get(paste0(NAM_list[i,2],tissue_list[7])),by="pan",all.x = T)
for (i in 1:26) ear_mat = merge(ear_mat,get(paste0(NAM_list[i,2],tissue_list[9])),by="pan",all.x = T)
for (i in 1:26) tassel_mat = merge(tassel_mat,get(paste0(NAM_list[i,2],tissue_list[10])),by="pan",all.x = T)
for (i in 1:26) root_mat = merge(root_mat,get(paste0(NAM_list[i,2],tissue_list[3])),by="pan",all.x = T)
for (i in c(1:26)[-c(4,5,6,8,9,16)]) embryo_mat = merge(embryo_mat,get(paste0(NAM_list[i,2],tissue_list[1])),by="pan",all.x = T)
for (i in c(1:26)[-22]) tip_mat = merge(tip_mat,get(paste0(NAM_list[i,2],tissue_list[8])),by="pan",all.x = T)
for (i in c(1:26)[-c(4,6,9)]) endosperm_mat = merge(endosperm_mat,get(paste0(NAM_list[i,2],tissue_list[2])),by="pan",all.x = T)

methylation_expression <- function(gene_list,swi1,swi2,methylation,expression,tissue,index = 2:27){
  expres_swi = merge(gene_list,expression,by="pan",all.x=T)[,index]
  methy_swi = merge(gene_list,methylation,by="pan",all.x=T)[,index]
  express1 <- express2 <- expres_swi
  express1[methy_swi != swi1] <- NA
  express2[methy_swi != swi2] <- NA
  express1$tpm = rowMeans(as.matrix(express1),na.rm = T)
  express2$tpm = rowMeans(as.matrix(express2),na.rm = T)
  df = data.frame(express1$tpm,express2$tpm)
  df = cbind(df,tissue)
 # df = drop_na(data.frame(df))
  names(df) <- c(swi1,swi2,"Tissue")
  return(df)
}

df_UM_gbM_common = data.frame(pan=Matrix_epiallele$pan[Matrix_epiallele$status=="UM_gbM"&Matrix_epiallele$copy=="1-1"])
df_UM_gbM <- data.frame(UM=numeric(),gbM=numeric(),tissue = character())
for (j in c(3:7,9,10))  df_UM_gbM <- rbind(df_UM_gbM,methylation_expression(df_UM_gbM_common,"UM","gbM",single_methylation_matrix,get(paste0(tissue_list[j],"_mat")),tissue_list[j]))
embryo_UM_gbM <- methylation_expression(df_UM_gbM_common,"UM","gbM",single_methylation_matrix[,-c(4,5,6,8,9,16)],embryo_mat,tissue_list[1],index = 2:21)
tip_UM_gbM <- methylation_expression(df_UM_gbM_common,"UM","gbM",single_methylation_matrix[,-c(22)],tip_mat,tissue_list[8],index = 2:26)
endosperm_UM_gbM <- methylation_expression(df_UM_gbM_common,"UM","gbM",single_methylation_matrix[,-c(4,6,9)],endosperm_mat,tissue_list[2],index = 2:24)
df_UM_gbM = rbind(df_UM_gbM,embryo_UM_gbM);df_UM_gbM = rbind(df_UM_gbM,tip_UM_gbM);df_UM_gbM = rbind(df_UM_gbM,endosperm_UM_gbM)

df_UM_gbM$difference = df_UM_gbM$gbM - df_UM_gbM$UM
df_UM_gbM = drop_na(df_UM_gbM)
df_UM_gbM$Tissue= factor(df_UM_gbM$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))
#Plot FIg4A
pdf("/Users/x/Desktop/Plot/Gene_Methylation/Fig4A_difference_violin_no_text.pdf",
    width = 4,height = 3)
ggplot(df_UM_gbM,aes(x=Tissue,y=difference)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.size = 0.1) +
  coord_cartesian(ylim = c(-15,15)) +
  theme_bw() + theme(text = element_blank())
dev.off()

df_UM_gbM_proportion = data.frame(Npan=df_UM_gbM %>% group_by(Tissue) %>% summarise(length(difference)),
                                  gbM_UM = df_UM_gbM %>% group_by(Tissue) %>% summarise(sum(difference>0)),
                                  UM_gbM = df_UM_gbM %>% group_by(Tissue) %>% summarise(sum(difference<0)),
                                  tie = df_UM_gbM %>% group_by(Tissue) %>% summarise(sum(difference==0)))[,-c(3,5,7)]
names(df_UM_gbM_proportion) = c("Tissue","Total","gbM","UM","equal")
df_UM_gbM_proportion = data.frame(Tissue=rep(df_UM_gbM_proportion$Tissue,2),
                                  Proportion = c(df_UM_gbM_proportion$gbM/df_UM_gbM_proportion$Total,df_UM_gbM_proportion$UM/df_UM_gbM_proportion$Total),
                                  Category = rep(c("gbM>UM","gbM<UM"),each=10))


#########other 2 controls: UM-teM, gbM-teM
df_UM_teM_common <- data.frame(pan=Matrix_epiallele$pan[Matrix_epiallele$status=="UM_teM"&Matrix_epiallele$copy=="1-1" ])
df_UM_teM <- data.frame(UM=numeric(),teM=numeric(),tissue = character())
for (j in c(3:7,9,10))  df_UM_teM <- rbind(df_UM_teM,methylation_expression(df_UM_teM_common,"UM","teM",single_methylation_matrix,get(paste0(tissue_list[j],"_mat")),tissue_list[j]))
embryo_UM_teM <- methylation_expression(df_UM_teM_common,"UM","teM",single_methylation_matrix[,-c(4,5,6,8,9,16)],embryo_mat,tissue_list[1],index = 2:21)
tip_UM_teM <- methylation_expression(df_UM_teM_common,"UM","teM",single_methylation_matrix[,-c(22)],tip_mat,tissue_list[8],index = 2:26)
endosperm_UM_teM <- methylation_expression(df_UM_teM_common,"UM","teM",single_methylation_matrix[,-c(4,6,9)],endosperm_mat,tissue_list[2],index = 2:24)
df_UM_teM = rbind(df_UM_teM,embryo_UM_teM);df_UM_teM = rbind(df_UM_teM,tip_UM_teM);df_UM_teM = rbind(df_UM_teM,endosperm_UM_teM)
df_UM_teM$difference = df_UM_teM$UM - df_UM_teM$teM
df_UM_teM = drop_na(df_UM_teM)
df_gbM_teM_common <- data.frame(pan=Matrix_epiallele$pan[Matrix_epiallele$status=="gbM_teM" & Matrix_epiallele$copy=="1-1"])
df_gbM_teM <- data.frame(gbM=numeric(),teM=numeric(),tissue = character())
for (j in c(3:7,9,10))  df_gbM_teM <- rbind(df_gbM_teM,methylation_expression(df_gbM_teM_common,"gbM","teM",single_methylation_matrix,get(paste0(tissue_list[j],"_mat")),tissue_list[j]))
embryo_gbM_teM <- methylation_expression(df_gbM_teM_common,"gbM","teM",single_methylation_matrix[,-c(4,5,6,8,9,16)],embryo_mat,tissue_list[1],index = 2:21)
tip_gbM_teM <- methylation_expression(df_gbM_teM_common,"gbM","teM",single_methylation_matrix[,-c(22)],tip_mat,tissue_list[8],index = 2:26)
endosperm_gbM_teM <- methylation_expression(df_gbM_teM_common,"gbM","teM",single_methylation_matrix[,-c(4,6,9)],endosperm_mat,tissue_list[2],index = 2:24)
df_gbM_teM = rbind(df_gbM_teM,embryo_gbM_teM);df_gbM_teM = rbind(df_gbM_teM,tip_gbM_teM);df_gbM_teM = rbind(df_gbM_teM,endosperm_gbM_teM)

df_gbM_teM_proportion = data.frame(Npan=df_gbM_teM %>% group_by(Tissue) %>% summarise(length(difference)),
                                  teM_UM = df_gbM_teM %>% group_by(Tissue) %>% summarise(sum(difference>0,na.rm = T)),
                                  UM_teM = df_gbM_teM %>% group_by(Tissue) %>% summarise(sum(difference<0,na.rm = T)) )[,-c(3,5)]
colnames(df_gbM_teM_proportion) <- c("Tissue","Total","UM","teM")
df_gbM_teM_proportion = data.frame(Tissue=rep(df_gbM_teM_proportion$Tissue,2),
                                  Proportion = c(df_gbM_teM_proportion$UM/df_gbM_teM_proportion$Total,df_gbM_teM_proportion$teM/df_gbM_teM_proportion$Total),
                                  Category = rep(c("gbM>teM","gbM<teM"),each=10))
df_gbM_teM_proportion$Tissue= factor(df_gbM_teM_proportion$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))



df_UM_teM_proportion = data.frame(Npan=df_UM_teM %>% group_by(Tissue) %>% summarise(length(difference)),
                                  gbM_UM = df_UM_teM %>% group_by(Tissue) %>% summarise(sum(difference>0)),
                                  UM_gbM = df_UM_teM %>% group_by(Tissue) %>% summarise(sum(difference<0)),
                                  tie = df_UM_teM %>% group_by(Tissue) %>% summarise(sum(difference==0)))[,-c(3,5,7)]
names(df_UM_teM_proportion) = c("Tissue","Total","UM","teM","equal")

df_UM_teM_proportion = data.frame(Tissue=rep(df_UM_teM_proportion$Tissue,2),
                                   Proportion = c(df_UM_teM_proportion$UM/df_UM_teM_proportion$Total,df_UM_teM_proportion$teM/df_UM_teM_proportion$Total),
                                   Category = rep(c("UM>teM","UM<teM"),each=10))
df_UM_teM_proportion$Tissue= factor(df_UM_teM_proportion$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))


df_UM_gbM %>% group_by(Tissue) %>% summarise(wilcox.test(UM,gbM)$p.value)
df_UM_gbM %>% group_by(Tissue) %>% summarise(binom.test(sum(difference>0),sum(difference>0,difference<0))$p.value)

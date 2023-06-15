library(tidyverse)
library(cowplot)

#Read in the core gene and epiallele data set for the subsetting
B73.core<- read.csv("/Users/x/Desktop/Data/core/B73.core",
                    sep = "\t",
                    header = F,
                    col.names =  c("chr","start","end","strand","gene","pan","copy","duplicate"))[,5:6]

B73.epiallele <- read.csv("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                          sep = "\t",
                          header = F,
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))
B73.core.epi <- merge(B73.core,B73.epiallele,by="gene",all.x = T)


Hufford_data = read.table("/Users/x/Desktop/Data/expression/TPM/B73.tpm.txt",
                          header = T)
Hufford_data$status <- "tissue"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])>=1)==10] <- "expressed"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])<1)==10] <- "silenced"
table(Hufford_data$status)
df_Hufford_data_core  = merge(B73.core.epi,Hufford_data,by="gene",all.x = T)
table(df_Hufford_data_core$epiallele,df_Hufford_data_core$status)
df_Hufford_data_core <- df_Hufford_data_core[df_Hufford_data_core$epiallele!="ambiguous",]
df_2E = data.frame(epiallele=rep(c("gbM","teM","UM"),3),
                   status = rep(c("constitutive","silent","tissue-specific"),each=3), 
                   count = c(table(df_Hufford_data_core$epiallele,df_Hufford_data_core$status)))

df_2E$epiallele = factor(df_2E$epiallele,levels = c("UM","gbM","teM"))
df_2E$status = factor(df_2E$status,levels = c("constitutive","silent","tissue-specific"))
F2E1=ggplot(df_2E[df_2E$epiallele=="UM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")

F2E2=ggplot(df_2E[df_2E$epiallele=="gbM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")
F2E3 = ggplot(df_2E[df_2E$epiallele=="teM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")
F2E = plot_grid(F2E1,F2E2,F2E3,nrow = 1)

df_TPM_core = data.frame(tissue=rep(names(df_Hufford_data_core)[6:15],each=dim(df_Hufford_data_core)[1]),
                         TPM = unlist(c(df_Hufford_data_core[,6:15])),
                         epiallele = rep(df_Hufford_data_core$epiallele,10))
df_TPM_core$epiallele = factor(df_TPM_core$epiallele,levels = c("UM","gbM","teM"))
df_TPM_core$tissue = factor(df_TPM_core$tissue,levels = names(df_Hufford_data_core)[c(11,12,13,9,8,10,15,14,6,7)])
F2F = ggplot(df_TPM_core,aes(x=tissue,y=log10(TPM+1),fill=epiallele)) +
  geom_boxplot(outlier.size =0.01)+
  theme_bw() + ylab("log10(tpm+1)") +
  coord_cartesian(ylim=c(0,5)) +
  stat_boxplot(geom ='errorbar',position=position_dodge())

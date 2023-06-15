library(tidyverse)
library(extrafont)
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

CDS_TE = read.table("/Users/x/Desktop/Data/annotation/Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS_TE_intersect.bed")[,c(7,11)]
names(CDS_TE)<- c("gene","length")
CDS_TE_merge = aggregate(CDS_TE$length,by=list(CDS_TE$gene),sum)
names(CDS_TE_merge)<- c("gene","length")
df_CDS_TE_merge = merge(B73.core.epi,CDS_TE_merge,by="gene",all.x = T)
df_CDS_TE_merge$length[is.na(df_CDS_TE_merge$length)] = 0
(freq = table(df_CDS_TE_merge$length>100,df_CDS_TE_merge$epiallele))
df_CDS_TE = data.frame(epiallele =c("gbM","teM","UM"),
                       proportion = c(402/sum(7138,402),249/sum(249,471),
                                      426/sum(426,12244)))
df_CDS_TE$epiallele = factor(df_CDS_TE$epiallele,levels = c("UM","gbM","teM"))
F2D=ggplot(df_CDS_TE,aes(x=epiallele,y=proportion)) + theme_bw() +
  geom_bar(stat = "identity") 

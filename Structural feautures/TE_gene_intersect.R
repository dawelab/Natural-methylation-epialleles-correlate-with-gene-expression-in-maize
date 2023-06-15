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
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))[,9:10]
B73.core.epi <- merge(B73.core,B73.epiallele,by="gene",all.x = T)

B73_TE = read.table("/Users/x/Desktop/Data/TE/B73_TE_intersect.bed")[,c(5,9)]
table(B73_TE$V9)
B73_TE = B73_TE[B73_TE$V9!="LTR",]
B73_TE = B73_TE[!duplicated(B73_TE),]
names(B73_TE) = c("gene","TE")
B73_TE_core = merge(B73.core.epi,B73_TE,by="gene")
B73_TE_core <- B73_TE_core[B73_TE_core$epiallele!="ambiguous",]
df_B73_TE_core = table(B73_TE_core$epiallele,B73_TE_core$TE)/c(7541,710,12666)
df_B73_TE_core = data.frame(count=c(df_B73_TE_core),
                            epiallele=rep(row.names(df_B73_TE_core),9),
                            superfamily = rep(colnames(df_B73_TE_core),each=3)  )
df_B73_TE_core$epiallele = factor(df_B73_TE_core$epiallele,
                                  levels = c("UM","gbM","teM"))
df_B73_TE_core$superfamily = factor(df_B73_TE_core$superfamily,
                                    levels = c("CACTA","Harbinger","hAT","helitron","Mariner",
                                               "Copia","Gypsy","LINE","Mutator"))
#Plotting
ggplot(df_B73_TE_core ,aes(x=superfamily,y=count,fill = epiallele)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  xlab("") +  ylab("") +
  theme() #text = element_blank()

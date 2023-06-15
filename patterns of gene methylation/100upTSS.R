setwd("/Users/x/Desktop/R/TSS")
#Read in the TSS data
B73_TSS_mCG = read.table("B73.canno.100upTSS.mCG",header = F,
                       col.names = c("gene","level","context"))
B73_TSS_mCHG = read.table("B73.canno.100upTSS.mCHG",header = F, 
                       col.names = c("gene","level","context"))
B73_TSS_mCHH = read.table("B73.canno.100upTSS.mCHH",header = F,
                       col.names = c("gene","level","context"))
#Read in the B73 core gene set for subsetting
B73.core<- read.csv("/Users/x/Desktop/Data/core/B73.core",
                    sep = "\t",
                    header = F,
                    col.names =  c("chr","start","end","strand","gene","pan","copy","duplicate"))[,5:6]

B73.epiallele <- read.csv("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                          sep = "\t",
                          header = F,
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))[,9:10]
B73.core.epi <- merge(B73.core,B73.epiallele,by="gene",all.x = T)

B73_TSS_mCG_core <- merge(B73.core.epi,B73_TSS_mCG,by="gene",all.x = T)
B73_TSS_mCHG_core <- merge(B73.core.epi,B73_TSS_mCHG,by="gene",all.x = T)
B73_TSS_mCHH_core <- merge(B73.core.epi,B73_TSS_mCHH,by="gene",all.x = T)
B73_TSS <- rbind(B73_TSS_mCG_core,B73_TSS_mCHG_core,B73_TSS_mCHH_core)

#Remove the ambiguous genes 
df_TSS <- B73_TSS[B73_TSS$epiallele!="ambiguous",]
df_TSS$epiallele <- factor(df_TSS$epiallele,levels = c("UM","gbM","teM"))
df_TSS$context<- factor(df_TSS$context,levels = c("CG","CHG","CHH"))

df_TSS %>% ggplot(aes(x=epiallele,y=level,fill=context)) +
  geom_boxplot(outlier.size=.01, width = 0.5,
               position=position_dodge(width = 1)) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width = 1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1)) +
  theme_bw()

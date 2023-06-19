#Read data, there are 10 (Hufford et al.) + 4 Tissues (Warman et al)
library("tidyverse")
Hufford_data = read.table("/Users/x/Desktop/Data/expression/TPM/B73.tpm.txt",
                          header = T)
Hufford_data$status <- "tissue"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])>=1)==10] <- "expressed"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])<1)==10] <- "silenced"
table(Hufford_data$status)
#Import and merge the core gene set of B73
B73.core<- read.csv("/Users/x/Desktop/Data/core_gene/B73.class.txt",
                    sep = ",",
                    header = T) %>% filter(class=="Core Gene")

B73.epiallele <- read.csv("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                          sep = "\t",
                          header = F,
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))[,9:10]
B73.core.epi <- merge(B73.core,B73.epiallele,by="gene",all.x = T)

B73_core_epiallele = merge(B73.core,B73.epiallele,by="gene",all.x = T)
B73_core_epiallele = B73_core_epiallele[B73_core_epiallele$epiallele!="ambiguous",]
temp = merge(Hufford_data,B73_core_epiallele,by="gene") 
table(temp$status,temp$epiallele)
df_Hufford = data.frame(tpm = unlist(temp[,5:14],use.names = F),
                        tissue = rep(names(temp)[5:14],each=length(temp$gene)),
                        epiallele = rep(temp$epiallele,10))


df_Hufford %>% filter(tissue=="ear") %>% group_by(epiallele) %>% summarise(mean(tpm))
df_Hufford$tpm = df_Hufford$tpm +1
#df_Hufford<- df_Hufford[df_Hufford$epiallele!="ambiguous",]
df_Hufford$epiallele <- factor(df_Hufford$epiallele,levels = c("UM","gbM","teM"))
df_Hufford$tissue <- factor(df_Hufford$tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))
#pdf("Fig2F_10tissue_notext.pdf",width =8 ,height = 4)
F2F = ggplot(df_Hufford,aes(x=tissue,y=log10(tpm),fill=epiallele)) +
  geom_boxplot(outlier.size =0.01)+
  theme_bw() + ylab("log10(tpm+1)") +
  coord_cartesian(ylim=c(0,5)) +
  stat_boxplot(geom ='errorbar',position=position_dodge())# +
  theme(text = element_blank())

#expression pattern

df_pattern = data.frame(methylation = rep(c("UM","gbM","teM"),3),
                        pattern = rep(c("expressed","silenced","tissue"),each=3),
                        count =c(unlist(table(temp$epiallele,temp$status)[c(3,1,2),],use.names = F)) )
df_pattern$methylation = factor(df_pattern$methylation,levels = c("UM","gbM","teM"))
df_patternUM = df_pattern[df_pattern$methylation=="UM",]
df_patternUM$persent = df_patternUM$count/12666
df_patterngbM = df_pattern[df_pattern$methylation=="gbM",]
df_patterngbM$persent = df_patterngbM$count/7541
df_patternteM = df_pattern[df_pattern$methylation=="teM",]
df_patternteM$persent = df_patternteM$count/710
UM_pattern = ggplot(df_patternUM,aes(x=methylation,y=count,fill=pattern)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(~methylation) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")

  
gbM_pattern = ggplot(df_patterngbM,aes(x=methylation,y=count,fill=pattern)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(~methylation) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")

teM_pattern = ggplot(df_patternteM,aes(x=methylation,y=count,fill=pattern)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(~methylation) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")

teM_pattern2 = ggplot(df_patternteM,aes(x=methylation,y=count,fill=pattern)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(~methylation) +
  theme_bw() + theme(axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")
pdf("Fig2_pattern.pdf",width =8 ,height = 4)
plot_grid(UM_pattern,gbM_pattern,teM_pattern,teM_pattern2,nrow = 1)
dev.off()

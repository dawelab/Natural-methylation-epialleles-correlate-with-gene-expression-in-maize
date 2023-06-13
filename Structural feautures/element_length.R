setwd("/Users/x/Desktop/R/element_length")
#Read in the input data for element length
B73_exon_intron_CDS_UTR <- read.csv("B73.intron.exon.CDS.UTR.txt",
                                    sep = "\t",
                                    header = F,
                                    col.names = c("element","length","gene"))


B73_TE <- read.csv("/Users/x/Desktop/R/element_length/Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.txt",
                   sep = "\t",
                   header = F,
                   col.names = c("element","length","gene"))

#For different element, perform different operation 
#Get the individual UTR length for different genes
UTR_length <- B73_exon_intron_CDS_UTR %>% filter(element=="five_prime_UTR" |
                                                 element=="three_prime_UTR")
#Add up the element length for exon and intron
B73_exon <- B73_exon_intron_CDS_UTR %>% filter(element == "exon") %>%
              group_by(gene) %>% summarise(sum(length)) %>% as.data.frame() 
exon_length <- data.frame(element = "exon",
                          length= B73_exon$`sum(length)`,
                            gene= B73_exon$gene)

B73_intron <- B73_exon_intron_CDS_UTR %>% filter(element == "intron") %>%
  group_by(gene) %>% summarise(sum(length)) %>% as.data.frame()
intron_length <- data.frame(element = "intron",
                          length= B73_intron$`sum(length)`,
                          gene= B73_intron$gene)

B73_CDS <- B73_exon_intron_CDS_UTR %>% filter(element == "CDS") %>%
  group_by(gene) %>% summarise(sum(length)) %>% as.data.frame()
CDS_length <- data.frame(element = "CDS",
                            length= B73_CDS$`sum(length)`,
                            gene= B73_CDS$gene)
#Average the TE length in gene unit
B73_TE_mean <- aggregate(length~gene,
                         data = B73_TE,
                         FUN = mean)
TE_length <- data.frame(element="TE",
                        length=B73_TE_mean$length,
                        gene = B73_TE_mean$gene)


df_length <- rbind(UTR_length,exon_length,intron_length,CDS_length,TE_length)

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

#Subset the core gene and add the epiallele column
TE_core = merge(TE_length,B73.core.epi,by="gene",all.y = T)
TE_core$element[is.na(TE_core$element)] <- "TE"
TE_core$length[is.na(TE_core$length)] <- 0

UTR_core = merge(UTR_length,B73.core.epi,by="gene",all.y = T)
UTR_core$element[is.na(UTR_core$element)] <- "UTR"
UTR_core$length[is.na(UTR_core$length)] <- 0

intron_core = merge(intron_length,B73.core.epi,by="gene",all.y = T)
intron_core$element[is.na(intron_core$element)] <- "intron"
intron_core$length[is.na(intron_core$length)] <- 0

exon_core = merge(exon_length,B73.core.epi,by="gene",all.y = T)
exon_core$element[is.na(exon_core$element)] <- "exon"
exon_core$length[is.na(exon_core$length)] <- 0

CDS_core = merge(CDS_length,B73.core.epi,by="gene",all.y = T)
CDS_core$element[is.na(CDS_core$element)] <- "CDS"
CDS_core$length[is.na(CDS_core$length)] <- 0

df_length = rbind(CDS_core,exon_core,intron_core,TE_core,UTR_core)

df_length$element[df_length$element =="five_prime_UTR"| df_length$element =="three_prime_UTR"] <- "UTR"
df_length <- df_length[df_length$epiallele!="ambiguous",]
df_length$epiallele = factor(df_length$epiallele,levels = c("UM","gbM","teM"))
df_length$element = factor(df_length$element,
                                     levels = c("UTR","TE","CDS","exon","intron"))
ggplot(df_length,aes(x=epiallele,y=length,fill=element)) +
  geom_boxplot(outlier.size=.01,
               width = 0.5,
               position=position_dodge(width = 0.7)) +
  theme_bw() + 
  coord_cartesian(ylim=c(0,7000)) +
  xlab("") +
  stat_boxplot(geom ='errorbar',
               width = 0.5,
               position=position_dodge(width = 0.7)) + 
  scale_y_continuous(breaks = c(0,2000,4000,6000,8000)) + ylab("")

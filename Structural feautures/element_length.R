library(tidyverse)
setwd("/Users/x/Desktop/R/element_length")
#Read in the core gene and epiallele data set for the subsetting
B73.core<- read.csv("/Users/x/Desktop/Data/core_gene/B73.class.txt",
                    sep = ",",
                    header = T) %>% filter(class=="Core Gene")

B73.epiallele <- read.csv("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                          sep = "\t",
                          header = F,
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))[,9:10]
B73.core.epi <- merge(B73.core,B73.epiallele,by="gene",all.x = T)

#Read in the input data for element length
B73_exon_intron_CDS_UTR <- read.csv("B73_CDS_exon_intron_UTR.txt",
                                    sep = "\t",
                                    header = F,
                                    col.names = c("element","length","gene"))

B73_exon_intron_CDS_UTR$element[B73_exon_intron_CDS_UTR$element %in% c("five_prime_UTR","three_prime_UTR")] <- "UTR"

B73_TE <- read.csv("/Users/x/Desktop/R/element_length/Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.bed",
                   sep = "\t",
                   header = F,
                   col.names = c("chrx","startx","endx","x1","x2","strand","gene","chry","starty","endy","length"))
B73_TE <- B73_TE[,c("gene","length")]
B73_TE <- data.frame(element = "TE",B73_TE)
B73_genetic_length <- rbind(B73_exon_intron_CDS_UTR,B73_TE)

#Add up the genetic elements in gene unit
length_tidy <- function(ele){
  element_length <- B73_genetic_length[B73_genetic_length$element==ele,]
  element_length <- aggregate(element_length$length,by=list(element_length$gene),sum)
  element_length <- B73_genetic_length %>% filter(element == ele) %>%
    group_by(gene) %>% summarise(sum(length)) %>% as.data.frame() %>% unique()
  element_length$element <- ele
  names(element_length)[2] <- "length"
  element_length_core <- merge(B73.core.epi,element_length,all.x = T)
  element_length_core$length[is.na(element_length_core$length)] <- 0
  element_length_core$element[is.na(element_length_core$element)] <- ele
  element_length_core <- unique(element_length_core)
  return(element_length_core)
}

df_length <- rbind(length_tidy("UTR"),
                   length_tidy("TE"),
                   length_tidy("intron"),
                   length_tidy("exon"),
                   length_tidy("CDS"))


#Check the generated tidy data frame matched with the core gene or not and is there any na
table(df_length$element)
sum(is.na(df_length))
#Remove the ambiguous gene
df_length <- df_length[df_length$epiallele!="ambiguous",]
#Reorder the epiallele 
df_length$epiallele = factor(df_length$epiallele,levels = c("UM","gbM","teM"))
#Reorder the genetic elements
df_length$element = factor(df_length$element,
                           levels = c("UTR","TE","intron","CDS","exon"))

##plotting using ggplot2
ggplot(df_length,aes(x=epiallele,y=length,fill=element)) +
  geom_boxplot(outlier.size=0,
               width = 0.5,
               position=position_dodge(width = 0.7)) +
  theme_bw() + 
  coord_cartesian(ylim=c(0,7000)) +
  xlab("") +
  stat_boxplot(geom ='errorbar',
               width = 0.5,
               position=position_dodge(width = 0.7)) + 
  scale_y_continuous(breaks = c(0,2000,4000,6000,8000,10000,12000)) + ylab("") +
  theme(text = element_blank())

#calculate the average cumulative genetic element length
df_length %>% filter(element == "UTR")  %>% group_by(epiallele) %>% summarise(mean(length))
df_length %>% filter(element == "TE")  %>% group_by(epiallele) %>% summarise(mean(length))
df_length %>% filter(element == "intron")  %>% group_by(epiallele) %>% summarise(mean(length))
df_length %>% filter(element == "exon")  %>% group_by(epiallele) %>% summarise(mean(length))
df_length %>% filter(element == "CDS")  %>% group_by(epiallele) %>% summarise(mean(length))

#Comparison the average lengths of the UM/gbM/teM gene sets
element <- c("UM","gbM","teM")
context <- c("UTR","CDS","exon","intron","TE")
compare <- combn(element,2)
p_value <- c()
gr_context <- c()
gr_element <- c()
for (i in 1:3) {
  for (j in 1:5){
    p1 = wilcox.test(df_length$length[df_length$epiallele==compare[1,i] & df_length$element== context[j]],
                     df_length$length[df_length$epiallele==compare[2,i] & df_length$element== context[j]])$p.value
    p_value <- c(p1,p_value)
    gr_context <- c(gr_context,paste0(compare[1,i],compare[2,i]))
    gr_element <- c(gr_element,context[j])
  }
}
cbind(p_value,gr_context,gr_element)

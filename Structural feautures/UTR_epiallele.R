library(tidyverse)
library(extrafont)
#Read in the core gene and epiallele data set for the subsetting
B73.core<- read.table("/Users/x/Desktop/Data/core_gene/B73.class.txt",
                      header = T, sep = ",")

B73.epiallele <- read.csv("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                          sep = "\t",
                          header = F,
                          col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele"))
B73.core.epi <- merge(B73.core[B73.core$class=="Core Gene",],
                       B73.epiallele,by="gene",
                       all.x = T
                       )
table(B73.core.epi$epiallele)
B73_UTR <- read.csv("/Users/x/Desktop/R/element_length/B73_CDS_exon_intron_UTR.txt",
                    sep = "\t", header = F)
names(B73_UTR) <- c("element","length","gene")
B73_UTR$element[B73_UTR$element=="five_prime_UTR"] <-"UTR"
B73_UTR$element[B73_UTR$element=="three_prime_UTR"] <-"UTR"

element_length_UTR <- B73_UTR %>% filter(element == "UTR") %>% group_by(gene) %>%
  summarise(sum(length)) %>% as.data.frame()
names(element_length_UTR) <- c("gene","length")

element_length_UTR <- merge(B73.core.epi,element_length_UTR,all.x = T)
element_length_UTR$length[is.na(element_length_UTR$length)] = 0

element_length_UTR <- element_length_UTR %>% select(c("gene","epiallele","length")) %>% unique
table(element_length_UTR$epiallele)

(UTR_count = table(element_length_UTR$length>0,element_length_UTR$epiallele))
df_UTR = data.frame(proportion = c(11940/(11940+76),7511/(7511+29),291/(419+291)),
                    methylation = c("UM","gbM","teM"))
df_UTR$methylation = factor(df_UTR$methylation,levels = c("UM","gbM","teM"))
(F2C=ggplot(df_UTR,aes(x=methylation,y = proportion)) +
  geom_bar(stat = "identity") +
  theme_bw() )

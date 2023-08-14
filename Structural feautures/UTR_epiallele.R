library(tidyverse)
library(extrafont)
#Read in the core gene and epiallele data set for the subsetting
B73.core = read.csv("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Data/B73.all.csv")
table(B73.core$epiallele)
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
df_UTR = data.frame(proportion = c(11940/(11940+726),7511/(7511+29),291/(419+291)),
                    methylation = c("UM","gbM","teM"))
df_UTR$methylation = factor(df_UTR$methylation,levels = c("UM","gbM","teM"))
(F2C=ggplot(df_UTR,aes(x=methylation,y = proportion)) +
    geom_bar(stat = "identity") +
    theme_bw() )

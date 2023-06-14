#Read in the B73 canonical annotation file with just exon
#In terminal, we make use of the following code to generate the input file:
#awk '$3=="exon"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 | cut -f1 -d ';' | sed 's/Parent=//g' | sed 's/_T00.//g' | awk '{print $9}' > B73.exon.count
library(tidyverse)
B73_anno <- read.table("/Users/x/Desktop/Data/annotation/B73.exon.count",
                       col.names = "exon")
B73_exon_count <- aggregate(B73_anno$exon, by = list(B73_anno$exon), FUN=length)
names(B73_exon_count) <- c("gene","count")
#Read in B73 core gene set
B73.core <- read.table("/Users/x/Desktop/Data/core/B73.core") %>% select(c("V5")) %>% unique()
names(B73.core) <- c("gene")
#Merge the exon count and the gene ID to subset the core gene set in B73
B73_exon_count_core <- merge(B73.core, B73_exon_count, by = "gene", all.x = T)
#Read in the methylatin status in B73
B73_methyl <- read.table("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt") %>%
  select(c("V9","V10")) %>% unique()
names(B73_methyl) <- c("gene","epiallele")
B73_core_epi_exon <- merge(B73_exon_count_core,B73_methyl,by="gene", all.x = T) %>% unique() %>% filter(epiallele!="ambiguous")
B73_core_epi_exon$epiallele = factor(B73_core_epi_exon$epiallele,
                                     levels = c("UM","gbM","teM")) 
ggplot(B73_core_epi_exon, aes(x=epiallele, y = count)) +
  geom_boxplot(outlier.size = 0.1)  +
  theme_bw() +
  theme(legend.position = "none")

B73_core_epi_exon %>% group_by(epiallele) %>% summarise(mean(count))

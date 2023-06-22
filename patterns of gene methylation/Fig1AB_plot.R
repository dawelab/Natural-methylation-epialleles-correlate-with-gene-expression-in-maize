library(ggExtra) 
library(tidyverse)
#library(ggbreak)
#library (plotrix)
setwd("/Users/x/Desktop/Plot/Genic methylation")
B73_all = read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/epiallele/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                     sep = "\t", col.names = c('chr','start','end','strand',
                                               'mCG','cCG','mCHG','cCHG',
                                               'gene','epiallele'))
table(B73_all$epiallele)
B73_class = read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/B73.class.txt",
                       sep = ",",header = T)
B73_pan = merge(B73_class,B73_all,by="gene",all.y = T)
B73_core = B73_pan[B73_pan$class=="Core Gene",]

table(B73_pan$class[B73_pan$cCG >=40 & B73_pan$cCHG >=40],B73_pan$epiallele[B73_pan$cCG >=40 & B73_pan$cCHG >=40])[,c(4,2,3,1)]
table(B73_pan$epiallele[B73_pan$cCG >=40 & B73_pan$cCHG >=40])[c(4,2,3,1)]

#Calculate the frequency of B73 gene methylation types
table(B73_all$epiallele[B73_all$cCG >=40 & B73_all$cCHG>=40])
table(B73_core$epiallele[B73_core$cCG >=40 & B73_core$cCHG>=40])

f1=ggplot(B73_all[B73_all$cCG >=40 & B73_all$cCHG>=40,],aes(x = mCG, y = mCHG, color = epiallele)) + 
  geom_point(size = .01) + 
  theme_bw()   +
  xlab("CG methylation") +
  ylab("CHG methylation") +
  theme(text = element_text(family ="Arial"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank())+ 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  scale_x_continuous(breaks = seq(0,1,0.2),limits = c(0,1)) 

(f1 =  ggMarginal(f1, type = "histogram", groupColour = TRUE, groupFill = TRUE,bins=41)) 

f2 = ggplot(B73_core[B73_core$cCG >=40 & B73_core$cCHG>=40,],aes(x = mCG, y = mCHG, color = epiallele)) + 
  geom_point(size = .01) + 
  theme_bw() + 
  xlab("CG methylation") +
  ylab("CHG methylation") +
  scale_colour_discrete(breaks=c("UM","gbM","teM","ambiguous"))  +
  theme(text = element_text(family ="Arial"),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank())+ scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  scale_x_continuous(breaks = seq(0,1,0.2),limits = c(0,1))  +
  guides(shape = guide_legend(override.aes = list(size = 0.05))) 
(f2 =  ggMarginal(f2, type = "histogram", groupColour = TRUE, groupFill = TRUE,bins=41) )

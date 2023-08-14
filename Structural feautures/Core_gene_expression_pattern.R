library(tidyverse)
library(cowplot)

#Read in the core gene and epiallele data set for the subsetting
B73_expr <- read.csv("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Data/B73.all.csv") %>% filter(class == "Core Gene" & epiallele != "ambiguous")
B73_expr[B73_expr$gene == "Zm00001eb007510",] 
B73_expr <- B73_expr[,c(10,14:23)]
B73_expr$status <- "tissue"
B73_expr$status[rowSums(as.matrix(B73_expr[,2:11])>=1)==10] <- "expressed"
B73_expr$status[rowSums(as.matrix(B73_expr[,2:11])<1)==10] <- "silenced"
#There is one gene that are duplicated in the provided pangene table, we subtract that 1 fromt the following table
table(B73_expr$epiallele,B73_expr$status)
df_2E = data.frame(epiallele=rep(c("gbM","teM","UM"),3),
                   status = rep(c("constitutive","silent","tissue-specific"),each=3), 
                   count = c(table(B73_expr$epiallele,B73_expr$status)))

df_2E$epiallele = factor(df_2E$epiallele,levels = c("UM","gbM","teM"))
df_2E$status = factor(df_2E$status,levels = c("constitutive","silent","tissue-specific"))
F2E1=ggplot(df_2E[df_2E$epiallele=="UM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")

F2E2=ggplot(df_2E[df_2E$epiallele=="gbM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")
F2E3 = ggplot(df_2E[df_2E$epiallele=="teM",],aes(x=epiallele,y=count,fill=status)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y")  + facet_grid(~epiallele)+
  facet_grid(~epiallele) +
  theme_bw() + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab("")
F2E = plot_grid(F2E1,F2E2,F2E3,nrow = 1)

df_TPM_core = data.frame(tissue=rep(names(B73_expr)[2:11],each=dim(B73_expr)[1]),
                         TPM = unlist(c(B73_expr[,2:11])),
                         epiallele = rep(B73_expr$epiallele,10))
df_TPM_core$epiallele = factor(df_TPM_core$epiallele,levels = c("UM","gbM","teM"))
df_TPM_core$tissue = factor(df_TPM_core$tissue,levels = names(B73_expr)[c(9,8,7,4,5,10,6,11,3,2)])
F2F = ggplot(df_TPM_core,aes(x=tissue,y=log10(TPM+1),fill=epiallele)) +
  geom_boxplot(outlier.size =0.01)+
  theme_bw() + ylab("log10(tpm+1)") +
  coord_cartesian(ylim=c(0,5)) +
  stat_boxplot(geom ='errorbar',position=position_dodge())

library(ggExtra) 
library(tidyverse)
library(cowplot)
library(ggbreak)
setwd("/Users/x/Desktop/Plot/Gene_Methylation")
###############################################################################
#                             Genic pattern(scatter plot)                     #
###############################################################################
B73_all = read.table("/Users/x/Desktop/Data/Figure/Fig1/A-B/Zm-B73-REFERENCE-NAM-5.0.gene.mtr.ID.type.txt",sep = "\t")
B73_core = read.table("/Users/x/Desktop/Data/Figure/Fig1/A-B/Zm-B73-REFERENCE-NAM-5.0.gene.mtr.ID.strand.core",sep = "\t")
names(B73_all) <- names(B73_core) <- c('chr','strand','start','end',
                                       'cCG','mCG','cCHG','mCHG','cCHH','mCHH',
                                       'geneID','Gene Methylation Category')
B73_all=as.data.frame(B73_all[B73_all$cCG>=40&B73_all$cCHG>=40,])
B73_core=as.data.frame(B73_core[B73_core$cCG>=40&B73_core$cCHG>=40,])
#Calculate the frequency of B73 gene methylation types
table(B73_all$`Gene Methylation Category`)
table(B73_core$`Gene Methylation Category`)

plot_genic_pattern <- function(data,label,title){
  p = ggplot(data,aes(x = mCG, y = mCHG, color = `Gene Methylation Category`)) + 
    geom_point(size = 0.4) + 
    theme_bw() + 
    xlim(0,1) + ylim(0,1) +
    xlab("CG methylation level") +
    ylab("CHG methylation level") +
    scale_colour_discrete(name="Gene methylation category",
                          breaks=c("UM","gbM","teM","ambiguous"),
                          labels=label) +
    ggtitle(title) +
    theme(text = element_text(size=20),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          legend.position = c(0.4,0.83),
          legend.key = element_rect(fill = "transparent"),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_blank())+ 
    guides(shape = guide_legend(override.aes = list(size = 0.05))) 
  ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE,bins=41) 
}
all_genic=plot_genic_pattern(B73_all,c("UM(14,356)","gbM(8,173)","teM(3,391)", "ambiguous(4,387)"),"B73 all genes")
core_genic = plot_genic_pattern(B73_core,c("UM(12,624)","gbM(7,570)","teM(719)", "ambiguous(3,626)"),"B73 core genes")
pdf("/Users/x/Desktop/Plot/Gene_Methylation/Fig1A.pdf",width = 20, height = 9)
plot_grid(all_genic,core_genic)
dev.off()
###############################################################################
#                             Meta gene                                       #
###############################################################################
UM_mfg = t(read.table("/Users/x/Desktop/Data/Figure/Fig1/C-E/Zm-B73-REFERENCE-NAM-5.0.mfg.UM")[8,])
UM_mfg = UM_mfg[UM_mfg !="total_ave_mC"]
gbM_mfg = t(read.table("/Users/x/Desktop/Data/Figure/Fig1/C-E/Zm-B73-REFERENCE-NAM-5.0.mfg.gbM")[8,])
gbM_mfg = gbM_mfg[gbM_mfg !="total_ave_mC"]
teM_mfg = t(read.table("/Users/x/Desktop/Data/Figure/Fig1/C-E/Zm-B73-REFERENCE-NAM-5.0.mfg.teM")[8,])
teM_mfg = teM_mfg[teM_mfg !="total_ave_mC"]
MFG  = data.frame(Methylation_level =  as.numeric(c(UM_mfg,gbM_mfg,teM_mfg)), 
                  Context = rep(c("CG","CHG","CHH"), each = 100),
                  position = rep(c(1:100),3),
                  Category = rep(c("UM","gbM","teM"),each=300))
MFG$Category = factor(MFG$Category,levels=c("UM","gbM","teM"))
MFG_plot =   ggplot(MFG, aes(x=position, y = Methylation_level, color = Context )) +
  geom_line() + theme_bw() + ylim(0,1) +  #   
  geom_vline(xintercept = 50, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 30, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 70, color = "black", linetype = "dashed") + 
  xlab("Position")   +
  scale_x_continuous(breaks=c(0,30,70,100),
                     labels=c("-3k","TSS","TTS","+3k"),
                     position = "bottom") + 
  ylab("Methylation Level") +  
  theme(text = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA)) +
  guides(shape = guide_legend(override.aes = list(size = 0.05))) +
  facet_grid(~Category)
pdf("/Users/x/Desktop/Plot/Gene_Methylation/Fig1B.pdf",width = 15, height = 5)
plot_grid(MFG_plot)
dev.off()
###############################################################################
#                             Upstream TSS                                    #
###############################################################################
gbM_CG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.gbM.CG.TSS.100")
gbM_CHG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.gbM.CHG.TSS.100")
gbM_CHH = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.gbM.CHH.TSS.100")
UM_CG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.UM.CG.TSS.100")
UM_CHG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.UM.CHG.TSS.100")
UM_CHH = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.UM.CHH.TSS.100")
teM_CG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.teM.CG.TSS.100")
teM_CHG = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.teM.CHG.TSS.100")
teM_CHH = read.table("/Users/x/Desktop/Data/TSS/Zm-B73-REFERENCE-NAM-5.0.mtr.teM.CHH.TSS.100")
df1 = data.frame(Methylation_level = c(UM_CG$V6,gbM_CG$V6,teM_CG$V6,UM_CHG$V6,gbM_CHG$V6,teM_CHG$V6,UM_CHH$V6,gbM_CHH$V6,teM_CHH$V6), 
                 class = rep(c(rep("UM",length(UM_CG$V6)),
                               rep("gbM",length(gbM_CG$V6)),
                               rep("teM",length(teM_CG$V6))),3),
                 context = c(rep("CG",length(UM_CG$V6)+length(gbM_CG$V6)+length(teM_CG$V6)),
                             rep("CHG",length(UM_CG$V6)+length(gbM_CG$V6)+length(teM_CG$V6)),
                             rep("CHH",length(UM_CG$V6)+length(gbM_CG$V6)+length(teM_CG$V6))),
                 location=rep("First 100 upstream TSS",61977))
df1$class = factor(df1$class, levels = c("UM","gbM","teM"))
pdf("/Users/x/Desktop/Plot/Gene_Methylation/Fig1c.pdf",width = 5, height = 5)
ggplot(df1,aes(x=class,y=Methylation_level,fill = context)) +
  geom_boxplot(outlier.size=.01) +
  stat_boxplot(geom ='errorbar') + 
  theme_bw() +
  ylab("Methylation Level") +
  xlab("") + theme(text = element_text(size=20))
dev.off()
UM_mfg = t(read.table("/Users/x/Desktop/Data/mfg/B73.UM.mfg")[8,])
UM_mfg = UM_mfg[UM_mfg !="total_ave_mC"]
gbM_mfg = t(read.table("/Users/x/Desktop/Data/mfg/B73.gbM.mfg")[8,])
gbM_mfg = gbM_mfg[gbM_mfg !="total_ave_mC"]
teM_mfg = t(read.table("/Users/x/Desktop/Data/mfg/B73.teM.mfg")[8,])
teM_mfg = teM_mfg[teM_mfg !="total_ave_mC"]
MFG  = data.frame(Methylation =  as.numeric(c(UM_mfg,gbM_mfg,teM_mfg)), 
                  Context = rep(c("CG","CHG","CHH"), each = 100),
                  position = rep(c(1:100),3),
                  Category = rep(c("UM","gbM","teM"),each=300))
MFG$Category = factor(MFG$Category,levels=c("UM","gbM","teM"))
MFG_UM_left = MFG[MFG$Category=="UM" & MFG$position <=50,]
fc1 =  ggplot(MFG_UM_left, aes(x=position, y = Methylation, color = Context )) +
  geom_line()  +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 30, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=30) + 
  ylab("") +  xlab("") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category) + theme(legend.position = "none", 
                                axis.text.y = element_blank(),
                                axis.line.y.left = element_blank())


MFG_UM_right = MFG[MFG$Category=="UM" & MFG$position>50,]
fc2 =  ggplot(MFG_UM_right, aes(x=position, y = Methylation, color = Context )) +
  geom_line() +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 70, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=c(70)) + 
  ylab("") +  xlab("")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category)  + 
  theme(legend.position = "none", 
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank() ) 
  
fc = plot_grid(fc1,fc2)

MFG_gbM_left = MFG[MFG$Category=="gbM" & MFG$position <=50,]
fd1 =  ggplot(MFG_gbM_left, aes(x=position, y = Methylation, color = Context )) +
  geom_line()  +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 30, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=30) + 
  ylab("") +  xlab("") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category) + theme(legend.position = "none", axis.text.y = element_blank())
                            
MFG_gbM_right = MFG[MFG$Category=="gbM" & MFG$position>50,]
fd2 =  ggplot(MFG_gbM_right, aes(x=position, y = Methylation, color = Context )) +
  geom_line() +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 70, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=c(70)) + 
  ylab("") +  xlab("")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category)  + theme(legend.position = "none", 
                                 axis.text.y  = element_blank(),
                                 axis.ticks.y = element_blank()) 

fd = plot_grid(fd1,fd2)

MFG_teM_left = MFG[MFG$Category=="teM" & MFG$position <=50,]
fe1 =  ggplot(MFG_teM_left, aes(x=position, y = Methylation, color = Context )) +
  geom_line()  +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 30, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=30) + 
  ylab("") +  xlab("") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category) + theme(legend.position = "none", axis.text.y = element_blank())


MFG_teM_right = MFG[MFG$Category=="teM" & MFG$position>50,]
fe2 =  ggplot(MFG_teM_right, aes(x=position, y = Methylation, color = Context )) +
  geom_line() +
  theme_bw( base_family = "mono")  +
  geom_vline(xintercept = 70, color = "black", linetype = "dashed") + 
  scale_x_continuous(breaks=c(70)) + 
  ylab("") +  xlab("")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))  +
  facet_grid(~Category)  + theme(legend.position = "none", 
                                 axis.text.y  = element_blank(),
                                 axis.ticks.y = element_blank()) 

fe = plot_grid(fe1,fe2)

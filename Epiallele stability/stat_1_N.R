library(readr)
epi_matrix_final <- read_rds("/Users/x/Desktop/epi_matrix_final.rds")
#view(epi_matrix_final)
count_epi <- function(epiallele){
  mat <- epi_matrix_final %>% filter(status == epiallele) %>% select(names(epi_matrix_final)[2:27])
  mat[mat==epiallele] <- NA
  epi_count <- c(unlist(mat,use.names = F))
  return(table(epi_count[!is.na(epi_count)]))
}
stat_2 <- rbind(count_epi("UM"),
      count_epi("gbM"),
      count_epi("teM"))[,c(4,2,3)]
rownames(stat_2) <- c("UM","gbM","teM")
print(stat_2)
stat_2/rowSums(stat_2)

#Read in the copy matrix to filter out all the genes with copy number of 4
copy_matrix <- read_rds("/Users/x/Desktop/copy_matrix.rds")
df_pan = data.frame(pan=epi_matrix_final$pan)
new_copy_matrix <- merge(df_pan,copy_matrix, by = "pan", all.x = T)[,2:27]
epi_matrix2 <- epi_matrix_final[2:27]
epi_matrix2[new_copy_matrix<4] <- NA
epi_matrix2$pan <- df_pan
epi_matrix2$status <- epi_matrix_final$status
stat_4 <-
rbind(
epi_matrix2[epi_matrix2$status=="UM",names(copy_matrix)[2:27]] %>% unlist(use.names = F) %>% table,
epi_matrix2[epi_matrix2$status=="gbM",names(copy_matrix)[2:27]] %>% unlist(use.names = F) %>% table,
epi_matrix2[epi_matrix2$status=="teM",names(copy_matrix)[2:27]] %>% unlist(use.names = F) %>% table)[,c(4,2,3)]
rownames(stat_4) <- c("UM","gbM","teM")
stat_4/rowSums(stat_4)
chisq.test(rbind(stat_2[1,],stat_4[1,]))
chisq.test(rbind(stat_2[2,],stat_4[2,]))
fisher.test(rbind(stat_2[3,],stat_4[3,]))
##Plotting
df_pattern = data.frame(single = rep(row.names(stat_2),3),
                        multi = rep(colnames(stat_2),each=3),
                        count = unlist(c(stat_2/rowSums2(stat_2)),use.names = F) )
df_pattern$single = factor(df_pattern$single,levels = c("UM","gbM","teM"))
df_pattern$multi = factor(df_pattern$multi,levels = c("teM","gbM","UM"))

df_patternUM = df_pattern[df_pattern$single=="UM",]
df_patterngbM = df_pattern[df_pattern$single=="gbM",]
df_patternteM = df_pattern[df_pattern$single=="teM",]
(UM_pattern = ggplot(df_patternUM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))

(gbM_pattern = ggplot(df_patterngbM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))

(teM_pattern = ggplot(df_patternteM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))


(teM_pattern2 = ggplot(df_patternteM,aes(x=single,y=count,fill=multi)) +
  geom_bar(stat = "identity",position = "stack", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
  facet_grid(~single) +
  theme_bw() + theme(axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()) +
  xlab("") +ylab(""))

plot_grid(UM_pattern,gbM_pattern,teM_pattern,teM_pattern2,nrow = 1)


###Plotting
df_pattern = data.frame(single = rep(row.names(stat_4),3),
                        multi = rep(colnames(stat_4),each=3),
                        count = unlist(c(stat_4/rowSums2(stat_4)),use.names = F) )
df_pattern$single = factor(df_pattern$single,levels = c("UM","gbM","teM"))
df_pattern$multi = factor(df_pattern$multi,levels = c("teM","gbM","UM"))

df_patternUM = df_pattern[df_pattern$single=="UM",]
df_patterngbM = df_pattern[df_pattern$single=="gbM",]
df_patternteM = df_pattern[df_pattern$single=="teM",]
(UM_pattern = ggplot(df_patternUM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))

(gbM_pattern = ggplot(df_patterngbM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))

(teM_pattern = ggplot(df_patternteM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))


(teM_pattern2 = ggplot(df_patternteM,aes(x=single,y=count,fill=multi)) +
    geom_bar(stat = "identity",position = "stack", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
    facet_grid(~single) +
    theme_bw() + theme(axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    xlab("") +ylab(""))

plot_grid(UM_pattern,gbM_pattern,teM_pattern,teM_pattern2,nrow = 1)

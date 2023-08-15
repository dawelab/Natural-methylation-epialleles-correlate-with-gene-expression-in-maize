library(matrixStats)
library(tidyverse)
##Read in the pangene matrix downloaded from cyverse
pan_coreID <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv") %>% 
  filter(class == "Core Gene") %>% select("Pan_gene_ID")
names(pan_coreID) <- "pan"

#Read in the combine data frame 
path <- "/Users/x/Desktop/Data/combine/"
files <- list.files(path)

for (i in 1:26) {
  assign(
    gsub(".csv","",files[i]),
    read.csv(paste0(path,files[i]))
  )
}
NAM <- gsub(".all.csv","",files)
#Generate the copy matrix for 1-N / 1-1 filtering
for (i in 1:26) {
  assign(
    paste0(NAM[i],".copy"),
    get(gsub(".csv","",files[i])) %>% select(c("pan","copy")) %>% unique()
  )
}
##
copy_matrix <- merge(pan_coreID,B73.copy,by="pan",all.x = T)
for (i in 2:26) {
  copy_matrix <- merge(copy_matrix,
                       get( paste0(NAM[i],".copy")),
                       by = "pan",
                       all.x = T)
}
names(copy_matrix)  <- c("pan",NAM)

#assign the pangene status for copy matrix
copy_matrix$min <- rowMins(as.matrix(copy_matrix[,2:27]),na.rm = T)
copy_matrix$max <- rowMaxs(as.matrix(copy_matrix[,2:27]),na.rm = T)
copy_matrix$status <- "unknow"
copy_matrix$status[copy_matrix$min==1 & copy_matrix$max == 1] <- "1_1"
copy_matrix$status[copy_matrix$min==1 & copy_matrix$max > 1] <- "1_N"
copy_matrix$status[copy_matrix$min > 1 & copy_matrix$max > 1] <- "N_N"
table(copy_matrix$status)

#Generate the CDS matrix for median calculation
for (i in 1:26) {
  assign(
    paste0(NAM[i],".CDS"),
    get(gsub(".csv","",files[i])) %>% select(c("pan","length")) %>% unique()
  )
}
##Aggregate the multi-copy gene into one cell
for (i in 1:26) {
  assign(
    paste0(NAM[i],".CDS.agg"),
    aggregate(length~pan,get(paste0(NAM[i],".CDS")),FUN = paste0 %>% unlist)
  )
}
##
CDS_matrix <- merge(pan_coreID,B73.CDS.agg,by="pan",all.x = T)
for (i in 2:26) {
  CDS_matrix <- merge(CDS_matrix,
                      get(paste0(NAM[i],".CDS.agg")),
                      by = "pan",
                      all.x = T)
}
names(CDS_matrix)  <- c("pan",NAM)
##Make the multicopy gene as NA
sum(CDS_matrix$pan!=copy_matrix$pan)
sum(names(CDS_matrix)[1:27]!=names(copy_matrix)[1:27])
CDS_matrix_ij <- CDS_matrix[,2:27]
copy_matrix_ij <- copy_matrix[,2:27]
CDS_matrix_ij[copy_matrix_ij>1] <- NA
CDS_matrix_ij_num <-  as.matrix(CDS_matrix_ij) %>% as.numeric() %>% matrix(ncol=26)

pan_median <- rowMedians(CDS_matrix_ij_num,na.rm = T)
CDS_matrix_ij_num = as.data.frame(CDS_matrix_ij_num)
CDS_matrix_ij_num$pan = CDS_matrix$pan
CDS_matrix_ij_num$median = pan_median
names(CDS_matrix_ij_num) <- c(names(CDS_matrix)[2:27],"pan","median")
CDS_matrix_ij_num <- CDS_matrix_ij_num[,c(names(CDS_matrix),"median")]

##CDS matrix 

#Filter the genes with large CDS length difference
for (i in 1:dim(CDS_matrix_ij_num)[1]) {
  for (j in 2:27) {
    if ( !is.na(CDS_matrix_ij_num[i,j]) &
         abs((CDS_matrix_ij_num[i,j]-CDS_matrix_ij_num[i,28])/CDS_matrix_ij_num[i,28]) > 0.1)
    {
      CDS_matrix_ij_num[i,j] <- NA
    }
  }
}


####mCG & mCHG matrice
#######
for (i in 1:26) {
  assign(
    paste0(NAM[i],".mCG"),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 & mCHG <= 0.05 & copy == 1 ) %>% 
      select(c("pan","mCG")) %>% unique()
  )
}

mCG_matrix <- merge(pan_coreID,B73.mCG,by="pan",all.x = T)
for (i in 2:26) {
  mCG_matrix <- merge(mCG_matrix,
                      get(paste0(NAM[i],".mCG")),
                      by = "pan",
                      all.x = T)
}
names(mCG_matrix)  <- c("pan",NAM)
##Make the multicopy gene as NA
sum(mCG_matrix$pan!=CDS_matrix$pan)
sum(names(mCG_matrix)[1:27]!=names(mCG_matrix)[1:27])
#Remove the genes with big CDS difference in mCG matrix
mCG_matrix[is.na(CDS_matrix_ij_num[,-28])] <- NA
mCG_matrix$min <- rowMins(as.matrix(mCG_matrix[,2:27]),na.rm = T)
mCG_matrix$max <- rowMaxs(as.matrix(mCG_matrix[,2:27]),na.rm = T)
mCG_matrix$dif <- mCG_matrix$max - mCG_matrix$min
mCG_matrix[mCG_matrix$dif>0.2,] %>% dim()

mCG_matrix_filter <- mCG_matrix[mCG_matrix$dif>0.2,1:27] 
mCG_pan <- data.frame(pan=mCG_matrix_filter$pan)


for (i in 1:26) {
  assign(
    paste0(NAM[i],".mCHG"),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 &  copy == 1 ) %>% 
      select(c("pan","mCHG")) %>% unique()
  )
}

mCHG_matrix <- merge(pan_coreID,B73.mCHG,by="pan",all.x = T)
for (i in 2:26) {
  mCHG_matrix <- merge(mCHG_matrix,
                       get(paste0(NAM[i],".mCHG")),
                       by = "pan",
                       all.x = T)
}
names(mCHG_matrix)  <- c("pan",NAM)
##Make the multicopy gene as NA
sum(mCHG_matrix$pan!=CDS_matrix$pan)
sum(names(mCHG_matrix)[1:27]!=names(mCG_matrix)[1:27])
#Remove the genes with big CDS difference in mCG matrix
mCHG_matrix[is.na(CDS_matrix_ij_num[,-28])] <- NA
mCHG_matrix$min <- rowMins(as.matrix(mCHG_matrix[,2:27]),na.rm = T)
mCHG_matrix$max <- rowMaxs(as.matrix(mCHG_matrix[,2:27]),na.rm = T)
mCHG_matrix$dif <- mCHG_matrix$max - mCHG_matrix$min
mCHG_matrix[mCHG_matrix$dif>0.2,] %>% dim()
mCHG_matrix_filter <- mCHG_matrix[mCHG_matrix$dif>0.2,1:27] 
mCHG_pan <- data.frame(pan=mCHG_matrix_filter$pan)

#Read in the expression data
tpm_matrix <- function(tissue,panlist){
  for (i in 1:26) {
    assign(
      paste0(NAM[i],".", tissue),
      get(gsub(".csv","",files[i])) %>% 
        filter(cCG >=40 & cCHG >= 40 & mCHG <= 0.05 & copy == 1 ) %>% 
        select(c("pan",tissue)) %>% unique()
    )
  }
  
  tissue_matrix <- merge(pan_coreID,
                         get(paste0(NAM[1],".", tissue)),
                         ,by="pan",all.x = T)
  for (i in 2:26) {
    tissue_matrix <- merge(tissue_matrix,
                           get(paste0(NAM[i],".", tissue)),
                           by = "pan",
                           all.x = T)
  }
  names(tissue_matrix)  <- c("pan",NAM)
  ##Make the multicopy gene as NA
  sum(tissue_matrix$pan!=CDS_matrix$pan)
  sum(names(tissue_matrix)[1:27]!=names(mCHG_matrix)[1:27])
  #Remove the genes with big CDS difference in mCG matrix
  tissue_matrix[is.na(CDS_matrix_ij_num[,-28])] <- NA
  tissue_matrix_filter <- merge(panlist,tissue_matrix,
                                by = "pan", all.x = T)
  return(tissue_matrix_filter[,-1])
}

mCHG_anther_matrix <- tpm_matrix("anther",mCHG_pan)
mCHG_root_matrix <- tpm_matrix("root",mCHG_pan)
mCHG_shoot_matrix <- tpm_matrix("shoot",mCHG_pan)
mCHG_base_matrix <- tpm_matrix("base",mCHG_pan)
mCHG_middle_matrix <- tpm_matrix("middle",mCHG_pan)
mCHG_tip_matrix <- tpm_matrix("tip",mCHG_pan)
mCHG_ear_matrix <- tpm_matrix("ear",mCHG_pan)
mCHG_tassel_matrix <- tpm_matrix("tassel",mCHG_pan)

mCG_anther_matrix <- tpm_matrix("anther",mCG_pan)
mCG_root_matrix <- tpm_matrix("root",mCG_pan)
mCG_shoot_matrix <- tpm_matrix("shoot",mCG_pan)
mCG_base_matrix <- tpm_matrix("base",mCG_pan)
mCG_middle_matrix <- tpm_matrix("middle",mCG_pan)
mCG_tip_matrix <- tpm_matrix("tip",mCG_pan)
mCG_ear_matrix <- tpm_matrix("ear",mCG_pan)
mCG_tassel_matrix <- tpm_matrix("tassel",mCG_pan)

stat_cal <- function(mat_tissue,CG_matrix){
  cor_list <- c()
  p_list <- c()
  for (i in c(1:dim(mat_tissue)[1])){
    x = as.numeric(as.matrix(unlist(mat_tissue[i,])))
    y = as.numeric(as.matrix(unlist(CG_matrix[i,])))
    
    if(sum(!is.na(x))>=3 & sd(x,na.rm = T)!=0 & sd(y,na.rm = T)!=0 ) {
      stat = cor.test(x,y)$estimate
      p_val = cor.test(x,y)$p.value
      cor_list <- c(cor_list,stat)
      p_list <- c(p_list,p_val)
    }
    else {
      stat = NA
      p_val = NA
      cor_list <- c(cor_list,stat)
      p_list <- c(p_list,p_val)
    }
  }
  df = data.frame(cor=cor_list,p=p_list)
  return(df)
}

mCHG_anther <- stat_cal(mCHG_anther_matrix,mCHG_matrix_filter[,2:27])
mCHG_root <- stat_cal(mCHG_root_matrix,mCHG_matrix_filter[,2:27])
mCHG_shoot <- stat_cal(mCHG_shoot_matrix,mCHG_matrix_filter[,2:27])
mCHG_base <- stat_cal(mCHG_base_matrix,mCHG_matrix_filter[,2:27])
mCHG_middle <- stat_cal(mCHG_middle_matrix,mCHG_matrix_filter[,2:27])
mCHG_tip <- stat_cal(mCHG_tip_matrix,mCHG_matrix_filter[,2:27])
mCHG_ear <- stat_cal(mCHG_ear_matrix,mCHG_matrix_filter[,2:27])
mCHG_tassel <- stat_cal(mCHG_tassel_matrix,mCHG_matrix_filter[,2:27])

mCG_anther <- stat_cal(mCG_anther_matrix,mCG_matrix_filter[,2:27])
mCG_root <- stat_cal(mCG_root_matrix,mCG_matrix_filter[,2:27])
mCG_shoot <- stat_cal(mCG_shoot_matrix,mCG_matrix_filter[,2:27])
mCG_base <- stat_cal(mCG_base_matrix,mCG_matrix_filter[,2:27])
mCG_middle <- stat_cal(mCG_middle_matrix,mCG_matrix_filter[,2:27])
mCG_tip <- stat_cal(mCG_tip_matrix,mCG_matrix_filter[,2:27])
mCG_ear <- stat_cal(mCG_ear_matrix,mCG_matrix_filter[,2:27])
mCG_tassel <- stat_cal(mCG_tassel_matrix,mCG_matrix_filter[,2:27])


###########embryo & endosperm do not have complete mRNA sequence data, make the matrice sperately
#Create embryo tpm matrix
#embryo lacking NAMs: CML52 CML228 CML277 CML247 CML333 M162W
tissue="embryo"
embryo_index <- c(1:26)[-which(NAM %in% c("CML52","CML228","CML277","CML247","CML333","M162W"))]
for (i in embryo_index) {
  assign(
    paste0(NAM[i],".", tissue),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 & mCHG <= 0.05 & copy == 1 ) %>% 
      select(c("pan",tissue)) %>% unique()
  )
}

embryo_matrix <- merge(pan_coreID,
                       get(paste0(NAM[1],".", tissue)),
                       by="pan",all.x = T)
for (i in embryo_index[-1]) {
  embryo_matrix <- merge(embryo_matrix,
                         get(paste0(NAM[i],".", tissue)),
                         by = "pan",
                         all.x = T)
}
names(embryo_matrix)  <- c("pan",NAM[embryo_index])
#Remove the genes with big CDS difference in mCG matrix
embryo_matrix[is.na(CDS_matrix_ij_num[,names(embryo_matrix)])] <- NA
mCHG_embryo_matrix_filter <- merge(mCHG_pan,embryo_matrix,
                                   by = "pan", all.x = T)
mCHG_embryo <- stat_cal(mCHG_embryo_matrix_filter[,-1],
                        mCHG_matrix_filter[,names(mCHG_embryo_matrix_filter[,-1])])

mCG_embryo_matrix_filter <- merge(mCG_pan,embryo_matrix,
                                  by = "pan", all.x = T)
mCG_embryo <- stat_cal(mCG_embryo_matrix_filter[,-1],
                       mCG_matrix_filter[,names(mCG_embryo_matrix_filter[,-1])])
#Create endosperm matrix
#CML52 CML228 CML277 does not have endosperm
endosperm_index <- c(1:26)[-which(NAM %in% c("CML52","CML228","CML277"))]
tissue="endosperm"
for (i in endosperm_index) {
  assign(
    paste0(NAM[i],".", tissue),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 & mCHG <= 0.05 & copy == 1 ) %>% 
      select(c("pan",tissue)) %>% unique()
  )
}

endosperm_matrix <- merge(pan_coreID,
                          get(paste0(NAM[1],".", tissue)),
                          by="pan",all.x = T)
for (i in endosperm_index[-1]) {
  endosperm_matrix <- merge(endosperm_matrix,
                            get(paste0(NAM[i],".", tissue)),
                            by = "pan",
                            all.x = T)
}
names(endosperm_matrix)  <- c("pan",NAM[endosperm_index])
#Remove the genes with big CDS difference in mCG matrix
endosperm_matrix[is.na(CDS_matrix_ij_num[,names(endosperm_matrix)])] <- NA
mCHG_endosperm_matrix_filter <- merge(mCHG_pan,endosperm_matrix,
                                      by = "pan", all.x = T)
mCG_endosperm_matrix_filter <- merge(mCG_pan,endosperm_matrix,
                                     by = "pan", all.x = T)
mCHG_endosperm <- stat_cal(mCHG_endosperm_matrix_filter[,-1],
                           mCHG_matrix_filter[,names(mCHG_endosperm_matrix_filter[,-1])])
mCG_endosperm <- stat_cal(mCG_endosperm_matrix_filter[,-1],
                          mCG_matrix_filter[,names(mCG_endosperm_matrix_filter[,-1])])

######Plotting
tissue <- c("anther","base","ear","middle","root","shoot","tassel","tip","endosperm","embryo")
(stat_c_mCG <-  data.frame(value = rbind(summary(mCG_anther$cor),
                                         summary(mCG_base$cor),
                                         summary(mCG_ear$cor),
                                         summary(mCG_middle$cor),
                                         summary(mCG_root$cor),
                                         summary(mCG_shoot$cor),
                                         summary(mCG_tassel$cor),
                                         summary(mCG_tip$cor),
                                         summary(mCG_endosperm$cor),
                                         summary(mCG_embryo$cor))[,3:4]  %>% as.list %>% unlist(use.names = F),
                           tissue= rep(tissue,2),
                           statistics = rep(c("Median","Mean"),each=10)))
stat_c_mCG$tissue = factor(stat_c$tissue,
                           levels =  c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo") )
stat_c_mCG$outlier <- "N"
stat_c_mCG$outlier[stat_c$tissue=="endosperm"] <- "Y" 
stat_c_mCG$outlier = factor(stat_c$outlier,
                            levels = c("Y","N"))


(stat_c_mCHG <-  data.frame(value = rbind(summary(mCHG_anther$cor),
                                          summary(mCHG_base$cor),
                                          summary(mCHG_ear$cor),
                                          summary(mCHG_middle$cor),
                                          summary(mCHG_root$cor),
                                          summary(mCHG_shoot$cor),
                                          summary(mCHG_tassel$cor),
                                          summary(mCHG_tip$cor),
                                          summary(mCHG_endosperm$cor),
                                          summary(mCHG_embryo$cor))[,3:4]  %>% as.list %>% unlist(use.names = F),
                            tissue= rep(tissue,2),
                            statistics = rep(c("Median","Mean"),each=10)))
stat_c_mCHG$tissue = factor(stat_c$tissue,
                            levels =  c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo") )
#######Supplemetal S4A
ear_cor_list = mCG_ear$cor
df_cor <- data.frame(cc=ear_cor_list)
ggplot(df_cor,aes(x=cc)) +
  geom_histogram(binwidth = 0.05,color="#88ada6") +
  theme_classic() +
  xlab("Correlation coefficient betweeen tpm & mCG") +
  ylab("Count") +
  ggtitle("                                              Ear")

##Supplemetal S4B
ggplot(stat_c,aes(x=tissue,y=value,shape=statistics,col=outlier)) +
  geom_point()+ 
  theme_classic()  +
  xlab("") +
  ylab("")


#######Supplemetal S4C
ear_cor_list = mCHG_ear$cor
df_cor <- data.frame(cc=ear_cor_list)
ggplot(df_cor,aes(x=cc)) +
  geom_histogram(binwidth = 0.05,color="#88ada6") +
  theme_classic() +
  xlab("Correlation coefficient betweeen tpm & mCG") +
  ylab("Count") +
  ggtitle("                              Ear")

##Supplemetal S4D
ggplot(stat_c_mCHG,aes(x=tissue,y=value,shape=statistics)) +
  geom_point(col="#00BFC4")+ 
  theme_classic()  +
  xlab("") +
  ylab("")

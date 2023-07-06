path <- "/Users/x/Desktop/Data/combine/"
files <- list.files(path)
NAM <- gsub(".all.csv","",files)
for (i in 1:26) {
  assign(
    gsub(".csv","",files[i]),
    read.csv(paste0(path,files[i]))
  )
}
#Read in the expression matrice data
tpm_matrix <- function(tissue,panlist){
  for (i in 1:26) {
    assign(
      paste0(NAM[i],".", tissue),
      get(gsub(".csv","",files[i])) %>% 
        filter(cCG >=40 & cCHG >= 40 & copy == 1 ) %>% 
        select(c("pan",tissue)) %>% unique()
    )
  }
  
  tissue_matrix <- merge(pan_coreID,
                         get(paste0(NAM[1],".", tissue)),
                         by="pan",all.x = T)
  for (i in 2:26) {
    tissue_matrix <- merge(tissue_matrix,
                           get(paste0(NAM[i],".", tissue)),
                           by = "pan",
                           all.x = T)
  }
  names(tissue_matrix)  <- c("pan",NAM)
  return(tissue_matrix)
}

anther_matrix <- tpm_matrix("anther")
root_matrix <- tpm_matrix("root")
shoot_matrix <- tpm_matrix("shoot")
base_matrix <- tpm_matrix("base")
middle_matrix <- tpm_matrix("middle")
tip_matrix <- tpm_matrix("tip")
ear_matrix <- tpm_matrix("ear")
tassel_matrix <- tpm_matrix("tassel")

tissue="embryo"
embryo_index <- c(1:26)[-which(NAM %in% c("CML52","CML228","CML277","CML247","CML333","M162W"))]
for (i in embryo_index) {
  assign(
    paste0(NAM[i],".", tissue),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40  & copy == 1 ) %>% 
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
endosperm_index <- c(1:26)[-which(NAM %in% c("CML52","CML228","CML277"))]
tissue="endosperm"
for (i in endosperm_index) {
  assign(
    paste0(NAM[i],".", tissue),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 & copy == 1 ) %>% 
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

#Read in the epiallele switch data
UM_gbM <-  read_rds(paste0(express_path,"epi_matrix_1_1_UM_gbM.rds"))
UM_teM <- read_rds(paste0(express_path,"epi_matrix_1_1_UM_teM.rds"))
gbM_teM <- read_rds(paste0(express_path,"epi_matrix_1_1_gbM_teM.rds"))

##Make the indication matrice for epiallele
#For UM_gbM switch matrix
UM_gbM_pan  <- data.frame(pan=UM_gbM$pan)
UM_gbM_UM <- matrix(0,nrow = dim(UM_gbM)[1],ncol = 26)
UM_gbM_gbM <- matrix(0,nrow = dim(UM_gbM)[1],ncol = 26)
UM_gbM_UM[UM_gbM[,2:27]=="UM"] <- 1
UM_gbM_gbM[UM_gbM[,2:27]=="gbM"] <- 1
#For UM_teM switch matrix
UM_teM_pan  <- data.frame(pan=UM_teM$pan)
UM_teM_UM <- matrix(0,nrow = dim(UM_teM)[1],ncol = 26)
UM_teM_teM <- matrix(0,nrow = dim(UM_teM)[1],ncol = 26)
UM_teM_UM[UM_teM[,2:27]=="UM"] <- 1
UM_teM_teM[UM_teM[,2:27]=="teM"] <- 1
#For gbM_teM switch matrix
gbM_teM_pan  <- data.frame(pan=gbM_teM$pan)
gbM_teM_gbM <- matrix(0,nrow = dim(gbM_teM)[1],ncol = 26)
gbM_teM_teM <- matrix(0,nrow = dim(gbM_teM)[1],ncol = 26)
gbM_teM_gbM[gbM_teM[,2:27]=="gbM"] <- 1
gbM_teM_teM[gbM_teM[,2:27]=="teM"] <- 1

##Subsetting tpm matrix corresponding to different switch matrix
tissue <- c("tip",  "middle", "base", "root", "shoot","ear",  "anther", "tassel","endosperm","embryo")

for (i in 1:10) {
  assign(
    paste0(tissue[i],"_matrix_UM_gbM"),
    merge(UM_gbM_pan,
          get(paste0(tissue[i],"_matrix")),
          all.x = T)[,-1]
  )
}

for (i in 1:10) {
  assign(
    paste0(tissue[i],"_matrix_UM_teM"),
    merge(UM_teM_pan,
          get(paste0(tissue[i],"_matrix")),
          all.x = T)[,-1]
  )
}

for (i in 1:10) {
  assign(
    paste0(tissue[i],"_matrix_gbM_teM"),
    merge(gbM_teM_pan,
          get(paste0(tissue[i],"_matrix")),
          all.x = T)[,-1]
  )
}

##Function for matrix operation
mat_operation <- function(indication1,indication2,tpm){
  indication1 <- as.matrix(indication1)
  indication2 <- as.matrix(indication2)
  
  identity <- matrix(1,nrow = dim(tpm)[1],ncol = dim(tpm)[2])
  tpm <- as.matrix(tpm)
  tpm[is.na(tpm)] <- 0
  diff1 <- diag( tpm %*% t(indication1))/diag(indication1 %*% t(indication1))
  diff2 <- diag( tpm %*% t(indication2))/diag(indication2 %*% t(indication2))
  df <- data.frame(diff1,diff2)
  return(df)
}

for (i in 1:8) {
  assign(
    paste0(tissue[i],".gbM_UM"),
    mat_operation(UM_gbM_gbM,UM_gbM_UM,
                  get(paste0(tissue[i],"_matrix_UM_gbM"))
                  
    ))
}

for (i in 1:8) {
  assign(
    paste0(tissue[i],".gbM_teM"),
    mat_operation(gbM_teM_gbM,gbM_teM_teM,
                  get(paste0(tissue[i],"_matrix_gbM_teM"))
                  
    ))
}


embryo.gbM_UM <- mat_operation(UM_gbM_gbM[,which(names(UM_gbM) %in% names(embryo_matrix_UM_gbM))-1],
                               UM_gbM_UM[,which(names(UM_gbM) %in% names(embryo_matrix_UM_gbM))-1],
                               embryo_matrix_UM_gbM)
endosperm.gbM_UM <- mat_operation(UM_gbM_gbM[,which(names(UM_gbM) %in% names(endosperm_matrix_UM_gbM))-1],
                                  UM_gbM_UM[,which(names(UM_gbM) %in% names(endosperm_matrix_UM_gbM))-1],
                                  endosperm_matrix_UM_gbM)

embryo.gbM_teM <- mat_operation(gbM_teM_gbM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                gbM_teM_teM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                embryo_matrix_gbM_teM)

embryo.gbM_teM <- mat_operation(gbM_teM_gbM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                gbM_teM_teM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                embryo_matrix_gbM_teM)


df_UM_gbM <- data_frame(tissue=character(),
                        diff = numeric())
for (i in 1:10) {
  df_UM_gbM <- rbind(df_UM_gbM,
                     data.frame(tissue=tissue[i],
                                diff = get(paste0(tissue[i],".gbM_UM"))))
}
df_UM_gbM$tissue= factor(df_UM_gbM$tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))

names(df_UM_gbM) <- c("tissue","gbM","UM")
df_UM_gbM$diff <- df_UM_gbM$gbM-df_UM_gbM$UM
df_UM_gbM %>% ggplot(aes(x=tissue,y=diff)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.size = 0.1) +
  coord_cartesian(ylim = c(-15,15)) +
  theme_bw() + theme(text = element_blank())

df_UM_gbM %>% group_by(tissue) %>% summarise(median(diff,na.rm = T)/mean(UM,na.rm = T))

df_UM_gbM_proportion = data.frame(Npan=df_UM_gbM %>% group_by(tissue) %>% summarise(length(diff)),
                                  gbM_UM = df_UM_gbM %>% group_by(tissue) %>% summarise(sum(diff>0,na.rm = T)),
                                  UM_gbM = df_UM_gbM %>% group_by(tissue) %>% summarise(sum(diff<0,na.rm = T)),
                                  tie = df_UM_gbM %>% group_by(tissue) %>% summarise(sum(diff==0,na.rm = T)))[,-c(3,5,7)]
names(df_UM_gbM_proportion) = c("Tissue","Total","gbM","UM","equal")
df_UM_gbM_proportion = data.frame(Tissue=rep(df_UM_gbM_proportion$Tissue,2),
                                  Proportion = c(df_UM_gbM_proportion$gbM/df_UM_gbM_proportion$Total,df_UM_gbM_proportion$UM/df_UM_gbM_proportion$Total),
                                  Category = rep(c("gbM>UM","gbM<UM"),each=10))

df_UM_gbM_proportion$Tissue= factor(df_UM_gbM_proportion$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))

df_UM_gbM_proportion %>% ggplot(aes(x=Tissue,y=Proportion,fill=Category)) +
  geom_bar(stat = "identity",
           position = "dodge") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

df_UM_gbM %>% group_by(tissue) %>%  summarise(
  binom.test(sum(diff>0,na.rm = T),
             sum(diff>0, diff<0,na.rm = T),
             alternative = "two.sided"
             )$p.value
)
df_UM_gbM %>% group_by(tissue) %>% summarise(wilcox.test(UM,gbM)$p.value)
#######
for (i in 1:8) {
  assign(
    paste0(tissue[i],".UM_teM"),
    mat_operation(UM_teM_UM,UM_teM_teM,
                  get(paste0(tissue[i],"_matrix_UM_teM"))
    ))
}
embryo.UM_teM <- mat_operation(UM_teM_UM[,which(names(UM_teM) %in% names(embryo_matrix_UM_teM))-1],
                               UM_teM_teM[,which(names(UM_teM) %in% names(embryo_matrix_UM_teM))-1],
                               embryo_matrix_UM_teM)

endosperm.UM_teM <- mat_operation(UM_teM_UM[,which(names(UM_teM) %in% names(endosperm_matrix_UM_teM))-1],
                                  UM_teM_teM[,which(names(UM_teM) %in% names(endosperm_matrix_UM_teM))-1],
                                  endosperm_matrix_UM_teM)

df_UM_teM <- data_frame(tissue=character(),
                        diff = numeric())
for (i in 1:10) {
  df_UM_teM <- rbind(df_UM_teM,
                     data.frame(tissue=tissue[i],
                                diff = get(paste0(tissue[i],".UM_teM"))))
}
df_UM_teM$tissue= factor(df_UM_teM$tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))
names(df_UM_teM) <- c("tissue","UM","teM")

df_UM_teM_proportion = data.frame(tissue,
                                  UM = (df_UM_teM %>% group_by(tissue) %>% summarise(sum(UM>teM, na.rm = T)))[,2],
                                  teM = (df_UM_teM %>% group_by(tissue) %>% summarise(sum(UM<teM, na.rm = T)))[,2],
                                  tie = (df_UM_teM %>% group_by(tissue) %>% summarise(sum(UM==teM, na.rm = T)))[,2]
                                  )
names(df_UM_teM_proportion) = c("Tissue","UM","teM","tie")
df_UM_teM_proportion$Total <- df_UM_teM_proportion$UM + df_UM_teM_proportion$teM + df_UM_teM_proportion$tie

df_UM_teM_proportion = data.frame(Tissue=rep(df_UM_teM_proportion$Tissue,2),
                                  Proportion = c(df_UM_teM_proportion$UM/df_UM_teM_proportion$Total,
                                                 df_UM_teM_proportion$teM/df_UM_teM_proportion$Total),
                                  Category = rep(c("UM>teM","UM<teM"),each=10))

df_UM_teM_proportion$Tissue= factor(df_UM_teM_proportion$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))

df_UM_teM_proportion %>% ggplot(aes(x=Tissue,y=Proportion,fill=Category)) +
  geom_bar(stat = "identity",
           position = "dodge") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


df_UM_gbM %>% group_by(tissue) %>%  summarise(binom.test(sum(diff>0,na.rm = T),sum(diff>0, diff<0,na.rm = T), alternative = "two.sided")$p.value)
df_UM_gbM %>% group_by(tissue) %>% summarise(wilcox.test(UM,gbM)$p.value)


df_UM_teM %>% group_by(tissue) %>%  summarise(binom.test(sum(UM>teM,na.rm = T),sum(UM>teM, UM<teM,na.rm = T), alternative = "two.sided")$p.value)
df_UM_teM %>% group_by(tissue) %>% summarise(wilcox.test(UM,teM)$p.value)

############
#######
for (i in 1:8) {
  assign(
    paste0(tissue[i],".gbM_teM"),
    mat_operation(gbM_teM_gbM,gbM_teM_teM,
                  get(paste0(tissue[i],"_matrix_gbM_teM"))
    ))
}
embryo.gbM_teM <- mat_operation(gbM_teM_gbM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                gbM_teM_teM[,which(names(gbM_teM) %in% names(embryo_matrix_gbM_teM))-1],
                                embryo_matrix_gbM_teM)

endosperm.gbM_teM <- mat_operation(gbM_teM_gbM[,which(names(gbM_teM) %in% names(endosperm_matrix_gbM_teM))-1],
                                   gbM_teM_teM[,which(names(gbM_teM) %in% names(endosperm_matrix_gbM_teM))-1],
                                  endosperm_matrix_gbM_teM)

df_gbM_teM <- data_frame(tissue=character(),
                        diff = numeric())
for (i in 1:10) {
  df_gbM_teM <- rbind(df_gbM_teM,
                     data.frame(tissue=tissue[i],
                                diff = get(paste0(tissue[i],".gbM_teM"))))
}
df_gbM_teM$tissue= factor(df_gbM_teM$tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))
names(df_gbM_teM) <- c("tissue","gbM","teM")

df_gbM_teM_proportion = data.frame(tissue,
                                  gbM = (df_gbM_teM %>% group_by(tissue) %>% summarise(sum(gbM>teM, na.rm = T)))[,2],
                                  teM = (df_gbM_teM %>% group_by(tissue) %>% summarise(sum(gbM<teM, na.rm = T)))[,2],
                                  tie = (df_gbM_teM %>% group_by(tissue) %>% summarise(sum(gbM==teM, na.rm = T)))[,2]
)
names(df_gbM_teM_proportion) = c("Tissue","gbM","teM","tie")
df_gbM_teM_proportion$Total <- df_gbM_teM_proportion$gbM + df_gbM_teM_proportion$teM + df_gbM_teM_proportion$tie

df_gbM_teM_proportion = data.frame(Tissue=rep(df_gbM_teM_proportion$Tissue,2),
                                  Proportion = c(df_gbM_teM_proportion$gbM/df_gbM_teM_proportion$Total,
                                                 df_gbM_teM_proportion$teM/df_gbM_teM_proportion$Total),
                                  Category = rep(c("gbM>teM","gbM<teM"),each=10))

df_gbM_teM_proportion$Tissue= factor(df_gbM_teM_proportion$Tissue,levels = c("tip","middle","base","root","shoot","ear","anther","tassel","endosperm","embryo"))

df_gbM_teM_proportion %>% ggplot(aes(x=Tissue,y=Proportion,fill=Category)) +
  geom_bar(stat = "identity",
           position = "dodge") + theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


df_UM_gbM %>% group_by(tissue) %>%  summarise(binom.test(sum(diff>0,na.rm = T),sum(diff>0, diff<0,na.rm = T), alternative = "two.sided")$p.value)
df_UM_gbM %>% group_by(tissue) %>% summarise(wilcox.test(UM,gbM)$p.value)


df_UM_teM %>% group_by(tissue) %>%  summarise(binom.test(sum(UM>teM,na.rm = T),sum(UM>teM, UM<teM,na.rm = T), alternative = "two.sided")$p.value)
df_UM_teM %>% group_by(tissue) %>% summarise(wilcox.test(UM,teM)$p.value)


file_names<- list.files("/Users/x/Desktop/Data/methylation/cgchgmtr")
path <- c("/Users/x/Desktop/Data/methylation/cgchgmtr/")
name <- gsub("Zm-","",file_names)
name <- gsub("-REFERENCE-NAM-1.0.1.canon.gene.gene.mtr.ID.type.txt","",name)
name <- gsub("-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt","",name)

for (i in 1:length(file_names)) {
  assign(paste0(name[i],"epiallele"), read.csv(paste0(path,file_names[i]), header = F,sep = "\t",
                        col.names = c("chr","start","end","strand","mCG","cCG",
                                      "mCHG","cCHG","gene","epiallele"))) }

file_names<- list.files("/Users/x/Desktop/Data/annotation/CDS")
for (i in 1:length(file_names)){
  assign(paste0(name[i],".CDS"), read.table(paste0("/Users/x/Desktop/Data/annotation/CDS/",file_names[i]), header = F,col.names = c("gene","length")))
  assign(paste0(name[i],".cumCDS"),aggregate(get(paste0(name[i],".CDS"))$length,by=list(get(paste0(name[i],".CDS"))$gene),sum))
}

file_names = list.files("/Users/x/Desktop/Data/pan_gene/version3")[-c(1,26)]

for (i in 1:26) {
  assign(paste0(name[i],"_class"),
         read.table(paste0("/Users/x/Desktop/Data/pan_gene/version3/",file_names[i]),
                    header = T,sep = ","))
}

names(B73.cumCDS) <- c("gene","length")
names(B97.cumCDS) <- c("gene","length")
names(CML103.cumCDS) <- c("gene","length")
names(CML228.cumCDS) <- c("gene","length")
names(CML247.cumCDS) <- c("gene","length")
names(CML277.cumCDS) <- c("gene","length")
names(CML322.cumCDS) <- c("gene","length")
names(CML333.cumCDS) <- c("gene","length")
names(CML52.cumCDS) <- c("gene","length")
names(CML69.cumCDS) <- c("gene","length")
names(HP301.cumCDS) <- c("gene","length")
names(Il14H.cumCDS) <- c("gene","length")
names(Ki11.cumCDS) <- c("gene","length")
names(Ki3.cumCDS) <- c("gene","length")
names(Ky21.cumCDS) <- c("gene","length")
names(M162W.cumCDS) <- c("gene","length")
names(M37W.cumCDS) <- c("gene","length")
names(Mo18W.cumCDS) <- c("gene","length")
names(Ms71.cumCDS) <- c("gene","length")
names(NC350.cumCDS) <- c("gene","length")
names(NC358.cumCDS) <- c("gene","length")
names(Oh43.cumCDS) <- c("gene","length")
names(Oh7B.cumCDS) <- c("gene","length")
names(P39.cumCDS) <- c("gene","length")
names(Tx303.cumCDS) <- c("gene","length")
names(Txi8.cumCDS) <- c("gene","length")
pan_cds = read.table("/Users/x/Desktop/Data/matrix/median_4_pan.txt",header = T)
names(pan_cds) <- c("pan","single_length")

for (i in 1:26){
  assign(paste0(name[i],"_pan_class"),
         merge(pan_cds,
         merge(get(paste0(name[i],"_class")),
         merge(get(paste0(name[i],"epiallele"))[,c("gene","epiallele")],
               get(paste0(name[i],".cumCDS")), by="gene"),by="gene"),
         by = "pan") %>% filter(length<=1.1*single_length|length>=0.9*single_length)
         %>% select(c("pan","epiallele"))) 
}

for (i in 1:26) {
  assign(paste0(name[i],"_epi"),
         aggregate(get(paste0(name[i],"_pan_class"))[,c("epiallele")],
         by = list(get(paste0(name[i],"_pan_class"))[,c("pan")]),
         paste
         ))
}

for (i in 1:26) {
  assign(paste0(name[i],"_copy"),
         aggregate(get(paste0(name[i],"_pan_class"))[,c("epiallele")],
                   by = list(get(paste0(name[i],"_pan_class"))[,c("pan")]),
                   length
         ))
}

pan_v3 = read.table("/Users/x/Desktop/Data/pan_gene/version3/pan_gene_matrix_v3_cyverse.csv",
                    sep = ",",header = T)
core_pan = data.frame(Group.1 =pan_v3$Pan_gene_ID[pan_v3$class=="Core Gene"])

Matrix_epiallele = merge(core_pan,B73_epi,by="Group.1",all.x = T)
for (i in 2:26) {
  Matrix_epiallele = merge(Matrix_epiallele,
                           get(paste0(name[i],"_epi")),
                           by="Group.1",all.x = T)
}

Matrix_epiallele[Matrix_epiallele=="ambiguous"] <- NA
names(Matrix_epiallele) <- c("pan",name)
#Matrix_epiallele = sapply(Matrix_epiallele, as.character)
#write.csv(Matrix_epiallele,
#          file = "/Users/x/Desktop/Data/matrix/Matrix_epiallele.csv",
#          quote = F,
#          row.names = F)



df_NAM_copy = merge(core_pan, B73_copy,by="Group.1",all.x = T)
for (i in 2:length(name)) df_NAM_copy = merge(df_NAM_copy,get(paste0(name[i],"_copy")),by="Group.1",all.x = T)
names(df_NAM_copy) <- c("pan",name)
sum(Matrix_epiallele$pan!=df_NAM_copy$pan)
sum(names(Matrix_epiallele)!=names(df_NAM_copy))
df_NAM_copy$Min = rowMins(as.matrix(df_NAM_copy[,2:27]),na.rm = T)
df_NAM_copy$Max = rowMaxs(as.matrix(df_NAM_copy[,2:27]),na.rm = T)
sum(df_NAM_copy$Min==1&df_NAM_copy$Max==1)
sum(df_NAM_copy$Min==1&df_NAM_copy$Max!=1)
sum(df_NAM_copy$Min!=1&df_NAM_copy$Max!=1)
Matrix_epiallele$copy = "unknown"
Matrix_epiallele$copy[df_NAM_copy$Min==1&df_NAM_copy$Max==1] = "1-1"
Matrix_epiallele$copy[df_NAM_copy$Min==1&df_NAM_copy$Max>1] = "1-N"

Matrix_epiallele$UM = rowSums(Matrix_epiallele[,2:27]=="UM")
Matrix_epiallele$gbM = rowSums(Matrix_epiallele[,2:27]=="gbM")
Matrix_epiallele$teM = rowSums(Matrix_epiallele[,2:27]=="teM")

Matrix_epiallele$status ="unknown"
Matrix_epiallele$status[Matrix_epiallele$UM>=2&Matrix_epiallele$gbM==0&Matrix_epiallele$teM==0] <- "UM"
Matrix_epiallele$status[Matrix_epiallele$UM==0&Matrix_epiallele$gbM>=2&Matrix_epiallele$teM==0] <- "gbM"
Matrix_epiallele$status[Matrix_epiallele$UM==0&Matrix_epiallele$gbM==0&Matrix_epiallele$teM>=2] <- "teM"
Matrix_epiallele$status[Matrix_epiallele$UM>=2&Matrix_epiallele$gbM>=2&Matrix_epiallele$teM==0] <- "UM_gbM"
Matrix_epiallele$status[Matrix_epiallele$UM>=2&Matrix_epiallele$gbM==0&Matrix_epiallele$teM>=2] <- "UM_teM"
Matrix_epiallele$status[Matrix_epiallele$UM==0&Matrix_epiallele$gbM>=2&Matrix_epiallele$teM>=2] <- "gbM_teM"
Matrix_epiallele$status[Matrix_epiallele$UM>=2&Matrix_epiallele$gbM>=2&Matrix_epiallele$teM>=2] <- "UM_gbM_teM"

table(Matrix_epiallele$status,Matrix_epiallele$copy)

Tandem_matrix <- Matrix_epiallele[Matrix_epiallele$copy=="1-N",]
Tandem_matrix <- Tandem_matrix[Tandem_matrix$status!="unknown",]
Tandem_matrix <- Tandem_matrix[Tandem_matrix$status!="gbM_teM",]
Tandem_matrix <- Tandem_matrix[Tandem_matrix$status!="UM_gbM",]
Tandem_matrix <- Tandem_matrix[Tandem_matrix$status!="UM_gbM_teM",]
Tandem_matrix <- Tandem_matrix[Tandem_matrix$status!="UM_teM",]
table(Tandem_matrix$status)
table(merge(Tandem_matrix,df_NAM_copy[,c("pan","Max")],by="pan")$Max)

single_multi_cal <- function(copy=4){
  multi_pan = data.frame(pan=df_NAM_copy$pan[df_NAM_copy$Max==copy & df_NAM_copy$Min==1])
  tandem_duplication_matrix = merge(multi_pan, Tandem_matrix, by ="pan",all.x = T)
  tandem_duplication_matrix_UM <- tandem_duplication_matrix[tandem_duplication_matrix$status=="UM",2:27];tandem_duplication_matrix_gbM <- tandem_duplication_matrix[tandem_duplication_matrix$status=="gbM",2:27];tandem_duplication_matrix_teM <- tandem_duplication_matrix[tandem_duplication_matrix$status=="teM",2:27]
  tandem_duplication_matrix_UM[tandem_duplication_matrix_UM=="UM"] <- NA
  tandem_duplication_matrix_gbM[tandem_duplication_matrix_gbM=="gbM"] <- NA
  tandem_duplication_matrix_teM[tandem_duplication_matrix_teM=="teM"] <- NA
  multi_methylation = rbind(table(c(unlist(tandem_duplication_matrix_UM))),table(c(unlist(tandem_duplication_matrix_gbM))),table(c(unlist(tandem_duplication_matrix_teM))))[,c(4,2,3)]
  row.names(multi_methylation) <- c("UM","gbM","teM")
  return(multi_methylation)
}
sin_mul_2 <- single_multi_cal(copy=2)
sin_mul_4 <- single_multi_cal(copy=4)
sin_mul_4[3,] <- c(0,0,144)
chi <- matrix(c(469,14863+105+256,144,144+1983+64),ncol = 2,nrow = 2,byrow = T)
chisq.test(chi,correct = F)

chi2 <- matrix(c(98,98+2625+469,58,76+622+58),ncol = 2,nrow = 2,byrow = T)
chisq.test(chi2,correct = F)
df_sin_mul_2 <- data.frame(proportion = unlist(c(sin_mul_2/rowSums(sin_mul_2))),
                           single = rep(c("UM","gbM","teM"),3),
                           multi = rep(c("UM","gbM","teM"),each=3))
df_sin_mul_4 <- data.frame(proportion = unlist(c(sin_mul_4/rowSums(sin_mul_4))),
                           single = rep(c("UM","gbM","teM"),3),
                           multi = rep(c("UM","gbM","teM"),each=3))
df_sin_mul_2$single <- factor(df_sin_mul_2$single,levels = c("UM","gbM","teM"))
df_sin_mul_4$single <- factor(df_sin_mul_4$single,levels = c("UM","gbM","teM"))
df_sin_mul_2$multi <- factor(df_sin_mul_2$multi,levels = c("UM","gbM","teM"))
df_sin_mul_4$multi <-factor(df_sin_mul_4$multi,levels = c("UM","gbM","teM"))

 

Matrix_epiallele$copy = "unknow"
Matrix_epiallele$copy[df_NAM_copy$Min==1&df_NAM_copy$Max>=2] = "1-N"
Matrix_epiallele$copy[df_NAM_copy$Min==1&df_NAM_copy$Max==1] = "1-1"
Matrix_epiallele$copy[df_NAM_copy$Min>=2&df_NAM_copy$Max>=2] = "N-N"
table(Matrix_epiallele$copy,Matrix_epiallele$status)
Fig3D = data.frame(proportion=c(4063/7061,1538/5953,290/379),
                   epiallele = c("UM","gbM","teM"))
#Check here
Fig3E = data.frame(proportion=c(7061/(7061+966+15+43),5953/(5953+32+966+15),379/(379+32+15+43)),
                   epiallele = c("UM","gbM","teM"))
Fig3D$epiallele = factor(Fig3D$epiallele,levels = c("UM","gbM","teM"))
Fig3E$epiallele = factor(Fig3E$epiallele,levels = c("UM","gbM","teM"))
table(Matrix_epiallele$status,Matrix_epiallele$copy)

##########################################################
#Calculate the epiallele number for B73 & 25 NAM founders#
##########################################################
name = name[-c(1:3)]
epiallele <- data.frame()
for (i in 1:26) {
  UM = table(get(paste0(name[i],"epiallele"))[,"epiallele"])["UM"]
  gbM = table(get(paste0(name[i],"epiallele"))[,"epiallele"])["gbM"]
  teM = table(get(paste0(name[i],"epiallele"))[,"epiallele"])["teM"]
  ambiguous = table(get(paste0(name[i],"epiallele"))[,"epiallele"])["ambiguous"]
  temp = data.frame(UM,gbM,teM,ambiguous)
  epiallele = rbind(epiallele,temp)
}
rownames(epiallele) <- name
epiallele$total <- rowSums(epiallele)
epiallele$pUM <- epiallele$UM/epiallele$total
epiallele$pgbM <- epiallele$gbM/epiallele$total
epiallele$pteM <- epiallele$teM/epiallele$total
epiallele$pambiguous <- epiallele$ambiguous/epiallele$total
epiallele[which.max(epiallele$pUM),]
epiallele[which.min(epiallele$pUM),]
epiallele[which.max(epiallele$pgbM),]
epiallele[which.min(epiallele$pgbM),]
epiallele[which.max(epiallele$pteM),]
epiallele[which.min(epiallele$pteM),]

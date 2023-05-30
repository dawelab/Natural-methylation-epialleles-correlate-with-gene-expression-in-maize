library(tidyverse)
warnings("off")
#Match the pangene and gene name to create the pan matrix for the CG/CHG methylation levels
#1.Read all core  gene files
coregene_path = "/Users/x/Desktop/Data/core/"
core_file <- list.files(coregene_path)
for (i in 1:length(core_file)) {
  assign(core_file[i],read_table(paste0(coregene_path,core_file[i]), col_names = 
                                   c("chr","start","end","strand","gene","pan","copy","duplicate"))[,5:8])
}

#2. Read all methylation level files
methy_path = "/Users/x/Desktop/Data/methylation/cgchgmtr/"
methy_file = list.files(methy_path)
methy_names = gsub("core","methy",core_file)
for (i in 1:length(methy_file)) {
  assign(methy_names[i],
         read.table(paste0(methy_path,methy_file[i]),
                    col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele")))
}
#3. Use merge to put together pangene names and gene names
pan_methy_names = gsub("core","pan_methy",core_file)
for (i in 1:length(pan_methy_names)) {
  assign(pan_methy_names[i],
         merge(get(core_file[i]),
               get(methy_names[i]),
               by = "gene")
  )
}
#4. Mask the ambiguous and teM gene mCG as NA and extract pan and mCG values to
#prepare the mCG pan gene matrix
mCG_list <- gsub("core","mCG",core_file)
for (i in 1:length(mCG_list)) {
  assign(mCG_list[i],
         select(filter(get(pan_methy_names[i]),
                       duplicate=="N"),
                c("pan","mCG","epiallele")
         )
  )
}
#5.mask the ambiguous and teM genes' mCG values
B73.mCG$mCG[B73.mCG$epiallele %in% c("ambiguous","teM")] <- NA
B97.mCG$mCG[B97.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML103.mCG$mCG[CML103.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML228.mCG$mCG[CML228.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML247.mCG$mCG[CML247.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML277.mCG$mCG[CML277.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML322.mCG$mCG[CML322.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML333.mCG$mCG[CML333.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML52.mCG$mCG[CML52.mCG$epiallele %in% c("ambiguous","teM")] <- NA
CML69.mCG$mCG[CML69.mCG$epiallele %in% c("ambiguous","teM")] <- NA
HP301.mCG$mCG[HP301.mCG$epiallele %in% c("ambiguous","teM")] <- NA
IL14H.mCG$mCG[IL14H.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Ki11.mCG$mCG[Ki11.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Ki3.mCG$mCG[Ki3.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Ky21.mCG$mCG[Ky21.mCG$epiallele %in% c("ambiguous","teM")] <- NA
M162W.mCG$mCG[M162W.mCG$epiallele %in% c("ambiguous","teM")] <- NA
M37W.mCG$mCG[M37W.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Mo18W.mCG$mCG[Mo18W.mCG$epiallele %in% c("ambiguous","teM")] <- NA
MS71.mCG$mCG[MS71.mCG$epiallele %in% c("ambiguous","teM")] <- NA
NC350.mCG$mCG[NC350.mCG$epiallele %in% c("ambiguous","teM")] <- NA
NC358.mCG$mCG[NC358.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Oh43.mCG$mCG[Oh43.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Oh7B.mCG$mCG[Oh7B.mCG$epiallele %in% c("ambiguous","teM")] <- NA
P39.mCG$mCG[P39.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Tx303.mCG$mCG[Tx303.mCG$epiallele %in% c("ambiguous","teM")] <- NA
Tzi8.mCG$mCG[Tzi8.mCG$epiallele %in% c("ambiguous","teM")] <- NA
#get rid of the epiallele column
for (i in 1:length(mCG_list)) {
  assign(mCG_list[i],
         select(get(pan_methy_names[i]),
                c("pan","mCG")
         )
  )
}
#aggregate the multiple copy genes
pre_methy_list <- gsub("methy","pre",methy_names)
for (i in 1:length(pre_methy_list)) {
  assign(pre_methy_list[i],
         aggregate(mCG~pan,
                   data = get(mCG_list[i]),
                   FUN = paste0%>%unlist)
  )
}

pan_mCG <- merge(get(pre_methy_list[1]),
                 get(pre_methy_list[2]),
                 by = "pan",
                 all = T)
for (i in 3:length(pre_methy_list)) {
  pan_mCG <- merge(pan_mCG,
                   get(pre_methy_list[i]),
                   by="pan",
                   all=T)
}
names(pan_mCG) <- c("pan",gsub(".core","",core_file))
#Filter out the mCG values with CDS length in big difference
NAM_CDS_matrix_int <- read.table("/Users/x/Desktop/Data/matrix/NAM_CDS_matrix.txt",                                     header = T)
mCG_pan_list = data.frame(pan=pan_mCG$pan)
CDS_matrix <- merge(mCG_pan_list,
                    NAM_CDS_matrix_int,
                    by = "pan",
                    all.x = T)
CDS_matrix<- CDS_matrix[,-28]
sum(names(pan_mCG) != names(CDS_matrix))
sum(pan_mCG$pan != CDS_matrix$pan)
pan_mCG[is.na(CDS_matrix)] <- NA
pan_mCG[pan_mCG=="NULL"] <- NA

pan_mCG2 <- as.matrix(pan_mCG[,2:27])  %>% as.numeric()
mat_mCG <- matrix(pan_mCG2,ncol = 26)

# Make expression data for all tissues
express <- "/Users/x/Desktop/Data/expression/TPM/"
express_file <- list.files(express)
express_list <- gsub(".txt","",express_file)
for (i in 1:26) {
  assign(express_list[i],
         read.table(paste0(express,express_file[i]),
                    header = T)
  )
}
#Create functions for generating TPM list for all tissues
#Read in tpm files for all NAM founders


TPM_matrix <- function(tissue){
  tissue <- as.character(tissue)
  #Create an two-column dataframe with gene name and tissue tpm values
  #and name as NAM_gene_list 
  #tissue <- "root"
  NAM_gene_list <- gsub("tpm","genelist",express_list)
  for (i in 1:26) {
    assign(NAM_gene_list[i],
           get(express_list[i]) %>% select(c("gene",tissue)))
  }
  
  #Match the panID with geneID
  for (i in 1:26) {
    assign(paste0(express_list[i],paste0(".",tissue)),
           select(merge(get(NAM_gene_list[i]),
                        get(core_file[i]),
                        by = "gene",
                        all.y = T),
                  c("pan",tissue))
    )  
  }
  #aggregate by pangene to combine multi-copy genes in one cell by panIDs
  for (i in 1:26) {
    assign(paste0(express_list[i],".pre"),
           aggregate(get(tissue)~pan,
                     data = get(paste0(express_list[i],".",tissue)),
                     FUN = paste0%>%unlist
           )
    )
  }
  #Create the TPM matrix
  tpm_matrix <- merge(get(paste0(express_list[1],".pre")),
                      get(paste0(express_list[2],".pre")),
                      by="pan",
                      all = T)
  for (i in 3:26) {
    tpm_matrix <- merge(tpm_matrix,
                        get(paste0(express_list[i],".pre")),
                        by="pan",
                        all = T
    )
  }
  names(tpm_matrix) <- c("pan",gsub(".core","",core_file))
  
  #Mask the genes wtih big CDS difference
  tpm <- merge(mCG_pan_list,
               tpm_matrix,
               by="pan",
               all.x = T)
  tpm[is.na(pan_mCG)] <- NA
  
  tpm2 <- as.matrix(tpm[,2:27]) %>% as.numeric()
  mat_tpm <- matrix(tpm2,ncol = 26)
  
  return(mat_tpm)
}
#generating the tpm matrix that mRNA sequence data is complete
root_matrix <- TPM_matrix("root")
middle_matrix <- TPM_matrix("middle")
shoot_matrix <-TPM_matrix("shoot")
anther_matrix <-TPM_matrix("anther")
base_matrix <- TPM_matrix("base")
tip_matrix <- TPM_matrix("tip")
ear_matrix <- TPM_matrix("ear")
tassel_matrix <- TPM_matrix("tassel")             

write.table(root_matrix,"/Users/x/Desktop/Data/matrix/root_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(middle_matrix,"/Users/x/Desktop/Data/matrix/middle_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(shoot_matrix,"/Users/x/Desktop/Data/matrix/shoot_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(anther_matrix,"/Users/x/Desktop/Data/matrix/anther_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(base_matrix,"/Users/x/Desktop/Data/matrix/base_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(tip_matrix,"/Users/x/Desktop/Data/matrix/tip_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(ear_matrix,"/Users/x/Desktop/Data/matrix/ear_matrix.txt",quote = F,
            row.names = F, col.names = F)
write.table(tassel_matrix,"/Users/x/Desktop/Data/matrix/tassel_matrix.txt",quote = F,
            row.names = F, col.names = F)

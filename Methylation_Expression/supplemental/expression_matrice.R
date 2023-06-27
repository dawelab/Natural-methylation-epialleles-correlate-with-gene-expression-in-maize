# Make expression data for all tissues
pan_matrix <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv")
(NAM <- names(pan_matrix)[4:29])

express <- "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Methylation_Expression/TPM/"
express_file <- paste0(NAM,".tpm.txt")
express_list <- gsub(".txt","",express_file)

for (i in 1:26) {
  assign(express_list[i],
         read.table(paste0(express,express_file[i]),
                    header = T)
  )
}
#Read in the core gene files
coregene_path = "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/"
core_file <- paste0(NAM,".class.txt")

for (i in 1:length(core_file)) {
  assign(core_file[i],read.table(paste0(coregene_path,core_file[i]), 
                                 header = T, sep = "," ) %>% filter(class == "Core Gene" & copy == 1))
}
#The core pangene ID for subsetting
core_pan_ID <- read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/core_panID.txt",
                          col.names = "pan")


##Read in the CDS matrix filtered out multi-copy genes and genes with big CDS difference
NAM_CDS_matrix_int <- read_rds("/Users/x/Desktop/NAM_CDS_matrix_int.rds") 
NAM_CDS_matrix_core <- merge(core_pan_ID,NAM_CDS_matrix_int,by="pan")[,1:27]

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
  tpm_matrix <- merge(core_panID,
                      get(paste0(express_list[1],".pre")),
                      by="pan",
                      all = T)
  for (i in 2:26) {
    tpm_matrix <- merge(tpm_matrix,
                        get(paste0(express_list[i],".pre")),
                        by="pan",
                        all = T
    )
  }
  names(tpm_matrix) <- c("pan",gsub(".core","",core_file))
  
  #Mask the genes wtih big CDS difference
  tpm <- merge(core_pan_ID,
               tpm_matrix,
               by="pan",
               all.x = T)
  tpm[is.na(NAM_CDS_matrix_core)] <- NA
  
  tpm2 <- as.matrix(tpm[,2:27]) %>% as.numeric()
  mat_tpm <- matrix(tpm2,ncol = 26) %>% as.data.frame()
  row.names(mat_tpm) <- core_pan_ID$pan
  colnames(mat_tpm) <- NAM
#  names(mat_tpm) <- NAM
#  mat_tpm$pan = core_pan_ID
 # mat_tpm = mat_tpm[,c("pan",NAM)]
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

write.csv(root_matrix,"/Users/x/Desktop/Data/matrix/root_matrix.csv", quote = F)
write.csv(middle_matrix,"/Users/x/Desktop/Data/matrix/middle_matrix.csv",quote = F)
write.csv(shoot_matrix,"/Users/x/Desktop/Data/matrix/shoot_matrix.csv",quote = F)
write.csv(anther_matrix,"/Users/x/Desktop/Data/matrix/anther_matrix.csv",quote = F)
write.csv(base_matrix,"/Users/x/Desktop/Data/matrix/base_matrix.csv",quote = F)
write.csv(tip_matrix,"/Users/x/Desktop/Data/matrix/tip_matrix.csv",quote = F)
write.csv(ear_matrix,"/Users/x/Desktop/Data/matrix/ear_matrix.csv",quote = F)
write.csv(tassel_matrix,"/Users/x/Desktop/Data/matrix/tassel_matrix.csv",quote = F)

###########embryo & endosperm do not have complete mRNA sequence data, make the matrice sperately
#Create embryo tpm matrix
#embryo lacking NAMs: CML52 CML228 CML277 CML247 CML333 M162W
NAM
tissue="embryo"
NAM_gene_list <- gsub("tpm","genelist",express_list)
for (i in c(1:3,5:13,15:16,20,22:26)) {
  assign(NAM_gene_list[i],
         get(express_list[i]) %>% select(c("gene",tissue)))
}

#Match the panID with geneID
for (i in c(1:3,5:13,15:16,20,22:26)) {
  assign(paste0(express_list[i],paste0(".",tissue)),
         select(merge(get(NAM_gene_list[i]),
                      get(core_file[i]),
                      by = "gene",
                      all.y = T),
                c("pan",tissue))
  )  
}
#aggregate by pangene to combine multi-copy genes in one cell by panIDs
for (i in  c(1:3,5:13,15:16,20,22:26)) {
  assign(paste0(express_list[i],".pre"),
         aggregate(get(tissue)~pan,
                   data = get(paste0(express_list[i],".",tissue)),
                   FUN = paste0%>%unlist
         )
  )
}
#Create the TPM matrix
embryo_tpm <- merge(core_pan_ID,
                    get(paste0(express_list[2],".pre")),
                    by="pan",
                    all.x = T)
for (i in  c(2:3,5:13,15:16,20,22:26)) {
  embryo_tpm <- merge(embryo_tpm,
                      get(paste0(express_list[i],".pre")),
                      by="pan",
                      all.x = T
  )
}
names(embryo_tpm) <- c("pan",NAM[c(1:3,5:13,15:16,20,22:26)])

#Mask the genes wtih big CDS difference
embryo_tpm[is.na(NAM_CDS_matrix_core[,names(embryo_tpm)])] <- NA
embryo_tpm2 <- as.matrix(embryo_tpm[,2:21]) %>% as.numeric()
embryo_tpm <- matrix(embryo_tpm2,ncol = 20)
colnames(embryo_tpm) <- NAM[c(1:3,5:13,15:16,20,22:26)]
rownames(embryo_tpm) <- core_pan_ID$pan
  
write.csv(embryo_tpm,"/Users/x/Desktop/Data/matrix/embryo_matrix.csv", quote = F)
#Create endosperm matrix
#CML52 CML228 CML277 does not have endosperm

tissue="endosperm"
for (i in c(1:13,15:16,18,20:26)) {
  assign(NAM_gene_list[i],
         get(express_list[i]) %>% select(c("gene",tissue)))
}

#Match the panID with geneID
for (i in c(1:13,15:16,18,20:26)) {
  assign(paste0(express_list[i],paste0(".",tissue)),
         select(merge(get(NAM_gene_list[i]),
                      get(core_file[i]),
                      by = "gene",
                      all.y = T),
                c("pan",tissue))
  )  
}
#aggregate by pangene to combine multi-copy genes in one cell by panIDs
for (i in c(1:13,15:16,18,20:26)) {
  assign(paste0(express_list[i],".pre"),
         aggregate(get(tissue)~pan,
                   data = get(paste0(express_list[i],".",tissue)),
                   FUN = paste0%>%unlist
         )
  )
}
#Create the TPM matrix
endosperm_tpm <- merge(core_pan_ID,
                       get(paste0(express_list[2],".pre")),
                       by="pan",
                       all.x = T)
for (i in c(2:13,15:16,18,20:26)) {
  endosperm_tpm <- merge(endosperm_tpm,
                         get(paste0(express_list[i],".pre")),
                         by="pan",
                         all.x = T
  )
}
names(endosperm_tpm) <- c("pan",NAM[c(1:13,15:16,18,20:26)])

#Mask the genes wtih big CDS difference
endosperm_tpm[is.na(NAM_CDS_matrix_core[,names(endosperm_tpm)])] <- NA
endosperm_tpm2 <- as.matrix(endosperm_tpm[,2:24]) %>% as.numeric()
endosperm_tpm <- matrix(endosperm_tpm2,ncol = 23)
colnames(endosperm_tpm) <- NAM[c(1:13,15:16,18,20:26)]
rownames(endosperm_tpm) <- core_pan_ID$pan
write.csv(endosperm_tpm, "/Users/x/Desktop/Data/matrix/endosperm_matrix.csv", quote = F)

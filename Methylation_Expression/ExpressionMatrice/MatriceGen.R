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
tissue_tpm_matrix <- function(tissue){
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
   return(tissue_matrix)
}
anther_matrix <- tissue_tpm_matrix("anther")
base_matrix <- tissue_tpm_matrix("base")
ear_matrix <- tissue_tpm_matrix("ear")
middle_matrix <- tissue_tpm_matrix("middle")
root_matrix <- tissue_tpm_matrix("root")
shoot_matrix <- tissue_tpm_matrix("shoot")
tassel_matrix <- tissue_tpm_matrix("tassel")
tip_matrix <- tissue_tpm_matrix("tip")


write.csv(anther_matrix, "/Users/x/Desktop/anther_matrix.csv",
          row.names = F,
          quote = F)
write.csv(base_matrix, "/Users/x/Desktop/base_matrix.csv",
          row.names = F,
          quote = F)
write.csv(ear_matrix, "/Users/x/Desktop/ear_matrix.csv",
          row.names = F,
          quote = F)
write.csv(middle_matrix, "/Users/x/Desktop/middle_matrix.csv",
          row.names = F,
          quote = F)
write.csv(root_matrix, "/Users/x/Desktop/root_matrix.csv",
          row.names = F,
          quote = F)
write.csv(shoot_matrix, "/Users/x/Desktop/shoot_matrix.csv",
          row.names = F,
          quote = F)
write.csv(tassel_matrix, "/Users/x/Desktop/tassel_matrix.csv",
          row.names = F,
          quote = F)
write.csv(tip_matrix, "/Users/x/Desktop/tip_matrix.csv",
          row.names = F,
          quote = F)


setwd("/Users/x/Desktop/Data/matrix/CDS/")
library(tidyverse)
library(matrixStats)
NAM_CDS_list = list.files()
NAM_CDS <- gsub("-REFERENCE-NAM-1.0.1.canon.cds_length.bed","",NAM_CDS_list)
NAM_CDS <- gsub("-REFERENCE-NAM-5.0.1.canon.cds_length.bed","",NAM_CDS)
NAM_CDS <- gsub("Zm-","",NAM_CDS)
NAM_CDS[26] <- "Tzi8"
#Import CDS length for each CDS element
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".CDS"),
         read.table(NAM_CDS_list[i],
                    col.names = c("gene","length"))
         )
}
#Add up the CDS length by genes
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".comCDS"),
         aggregate(length~gene,
                   data = get(paste0(NAM_CDS[i],".CDS")),
                   FUN = sum
                   )
        )
}
#Match the pan and gene ID 
coregene_path = "/Users/x/Desktop/Data/core_gene/"
core_file <- list.files(coregene_path)
for (i in 1:length(core_file)) {
  assign(core_file[i],read.table(paste0(coregene_path,core_file[i]), 
                                 sep = ",", header = T) %>%
           filter(class == "Core Gene")  ) 
}
#cbind(core_file,NAM_CDS_list)
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".preCDS"),
         merge(get(core_file[i]),
               get(paste0(NAM_CDS[i],".comCDS")),
               by="gene") %>% select(c("length","pan"))
        )
}
#aggregate the CDS matrix by pangenes
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".agrCDS"),
         aggregate(length~pan,
                   data = get(paste0(NAM_CDS[i],".preCDS")),
                   FUN = paste0%>%unlist
                   )
  )
}
#Create the  CDS matrix
NAM_CDS_matrix <- merge(get(paste0(NAM_CDS[1],".agrCDS")),
                        get(paste0(NAM_CDS[2],".agrCDS")),
                        by = "pan",
                        all = T)
for (i in 3:length(NAM_CDS)) {
  NAM_CDS_matrix <- merge(NAM_CDS_matrix,
                          get(paste0(NAM_CDS[i],".agrCDS")),
                          by = "pan",
                          all = T
  )
}
names(NAM_CDS_matrix) <- c("pan",gsub(".class.txt","",core_file))
#Create the copy matrix
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".copy"),
         aggregate(length~pan,
                   data = get(paste0(NAM_CDS[i],".preCDS")),
                   FUN = length
         )
  )
}

NAM_copy_matrix <- merge(get(paste0(NAM_CDS[1],".copy")),
                         get(paste0(NAM_CDS[2],".copy")),
                         by="pan",
                         all=T) 
for (i in 3:26) {
  NAM_copy_matrix <- merge(NAM_copy_matrix,
                           get(paste0(NAM_CDS[i],".copy")),
                           by="pan",
                           all=T)
}
names(NAM_copy_matrix) <- c("pan",gsub(".class.txt","",core_file))

#Compare index mCG matrix and copy number matrix
sum(names(NAM_copy_matrix) != names(NAM_CDS_matrix))
sum(NAM_copy_matrix$pan!=NAM_CDS_matrix$pan)
#Transform into double type
NAM_CDS_matrix[NAM_CDS_matrix=="NULL"] <- NA
NAM_CDS_matrix_int = NAM_CDS_matrix[,2:27]
NAM_copy_matrix_int = NAM_copy_matrix[,2:27]
NAM_CDS_matrix_int[NAM_copy_matrix_int>=2] <- NA
NAM_CDS_matrix_int <- as.matrix(NAM_CDS_matrix_int)
NAM_CDS_matrix_int <- as.numeric(NAM_CDS_matrix_int)
NAM_CDS_matrix_int <- matrix(NAM_CDS_matrix_int,
                             nrow = dim(NAM_copy_matrix_int)[1],
                             ncol = 26)
pan_median <- rowMedians(NAM_CDS_matrix_int,na.rm = T)
NAM_CDS_matrix_int = as.data.frame(NAM_CDS_matrix_int)
NAM_CDS_matrix_int$pan = NAM_CDS_matrix$pan
NAM_CDS_matrix_int$median = pan_median
names(NAM_CDS_matrix_int) <- c(names(NAM_CDS_matrix)[2:27],"pan","median")
NAM_CDS_matrix_int <- NAM_CDS_matrix_int[,c(names(NAM_CDS_matrix),"median")]

#Filter the genes with large CDS length difference
for (i in 1:dim(NAM_CDS_matrix_int)[1]) {
  for (j in 2:27) {
    if ( !is.na(NAM_CDS_matrix_int[i,j]) &
      abs((NAM_CDS_matrix_int[i,j]-NAM_CDS_matrix_int[i,28])/NAM_CDS_matrix_int[i,28]) > 0.1)
    {
      NAM_CDS_matrix_int[i,j] <- NA
    }
  }
}
write_delim(NAM_CDS_matrix_int,"/Users/x/Desktop/NAM_CDS_matrix.txt")

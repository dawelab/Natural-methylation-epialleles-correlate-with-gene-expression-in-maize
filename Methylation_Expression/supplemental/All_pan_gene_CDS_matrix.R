#In the canonical gene annotation files, use awk and (end - start + 1) to 
#calculate the CDS length for each element with gene ID
library(tidyverse)
library(matrixStats)
pan_matrix <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv")
(NAM <- names(pan_matrix)[4:29])
github_path <- "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Methylation_Expression/CDS/"
files <- c(paste0("Zm-",NAM[1],"-REFERENCE-NAM-5.0.1.canon.cds_length.bed"),
           paste0("Zm-",NAM[2:26],"-REFERENCE-NAM-1.0.1.canon.cds_length.bed"))

#Import CDS length for each CDS element
for (i in 1:length(NAM)) {
  assign(paste0(NAM[i],".CDS"),
         read.table(paste0(github_path,files[i]),
                    col.names = c("gene","length"))
  )
}
#Add up the CDS length by genes
for (i in 1:length(NAM)) {
  assign(paste0(NAM[i],".comCDS"),
         aggregate(length~gene,
                   data = get(paste0(NAM[i],".CDS")),
                   FUN = sum
         )
  )
}
#Match the pan and gene ID 
coregene_path = "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/"
core_file <- paste0(NAM,".class.txt")
for (i in 1:length(core_file)) {
  assign(paste0(NAM[i],".core"),
         read.table(paste0(coregene_path,core_file[i]), 
                                 sep = ",", header = T)  ) 
}
#cbind(core_file,NAM)
for (i in 1:length(NAM)) {
  assign(paste0(NAM[i],".preCDS"),
         merge(get(paste0(NAM[i],".core")),
               get(paste0(NAM[i],".comCDS")),
               by="gene") %>% select(c("length","pan"))
  )
}
#aggregate the CDS matrix by pangenes
for (i in 1:length(NAM)) {
  assign(paste0(NAM[i],".agrCDS"),
         aggregate(length~pan,
                   data = get(paste0(NAM[i],".preCDS")),
                   FUN = paste0%>%unlist
         )
  )
}
#Create the  CDS matrix
NAM_CDS_matrix <- merge(get(paste0(NAM[1],".agrCDS")),
                        get(paste0(NAM[2],".agrCDS")),
                        by = "pan",
                        all = T)
for (i in 3:length(NAM)) {
  NAM_CDS_matrix <- merge(NAM_CDS_matrix,
                          get(paste0(NAM[i],".agrCDS")),
                          by = "pan",
                          all = T
  )
}
names(NAM_CDS_matrix) <- c("pan",NAM)
#Create the copy matrix to mask all genes with copy number that is not 1
for (i in 1:length(NAM)) {
  assign(paste0(NAM[i],".copy"),
         aggregate(length~pan,
                   data = get(paste0(NAM[i],".preCDS")),
                   FUN = length
         )
  )
}

NAM_copy_matrix <- merge(get(paste0(NAM[1],".copy")),
                         get(paste0(NAM[2],".copy")),
                         by="pan",
                         all=T) 
for (i in 3:26) {
  NAM_copy_matrix <- merge(NAM_copy_matrix,
                           get(paste0(NAM[i],".copy")),
                           by="pan",
                           all=T)
}
names(NAM_copy_matrix) <- c("pan",NAM)

#Compare index CDS matrix and copy number matrix
sum(names(NAM_copy_matrix) != names(NAM_CDS_matrix))
sum(NAM_copy_matrix$pan!=NAM_CDS_matrix$pan)
#Transform into double type
#Subset the part without pangene information to do matrix replacement
NAM_CDS_matrix[NAM_CDS_matrix=="NULL"] <- NA

NAM_CDS_matrix_int = NAM_CDS_matrix[,2:27]
NAM_copy_matrix_int = NAM_copy_matrix[,2:27]

#replace the CDS matrix cell with copy number greater than 2
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

##Write out the CDS length matrix, here the Mulitple copy genes and genes with big
#difference in cumulatve CDS length is replaced with NA
write_rds(NAM_CDS_matrix_int,"/Users/x/Desktop/NAM_CDS_matrix_int.rds")

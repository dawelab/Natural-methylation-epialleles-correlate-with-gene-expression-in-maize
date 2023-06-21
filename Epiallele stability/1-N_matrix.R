library(tidyverse)
library(matrixStats)
#perform the filtering by cumulative CDS lengths
NAM_CDS_list = list.files("/Users/x/Desktop/Data/matrix/CDS")
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
  assign(gsub(".txt","",core_file[i]),read.table(paste0(coregene_path,core_file[i]),
                                                 sep = ",", header = T) %>%
           filter(class == "Core Gene")
  )
}
#cbind(core_file,NAM_CDS_list)
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".preCDS"),
         merge(get(core_file[i]),
               get(paste0(NAM_CDS[i],".comCDS")),
               by="gene", all.x = T) %>% select(c("pan","gene","length")) 
       
  )
}
#Read in the CDS matrix to extract the pangene and its median
CDS_matrix <- read.table("/Users/x/Desktop/NAM_CDS_matrix.txt",
                         header = T)[,c("pan","median")]
#Combine the median of CDS length together for filtering genes with 10% difference in CDS length
CDS_name <- gsub(".class.txt",".panCDS",core_file)
for (i in 1:26)
  assign(CDS_name[i],
         merge(CDS_matrix,get(paste0(NAM_CDS[i],".preCDS")),by="pan", all.x = T) %>% 
           filter(abs(length-median)/median <= 0.1))
}
#Subsetting the epialleles for filtered gene IDs
#Read in the methylation files
epi_path <- "/Users/x/Desktop/Data/methylation/cgchgmtr/"
epi_files <- list.files(epi_path)
#cbind(epi_files,core_file)
for (i in 1:26) {
  assign(gsub(".class.txt",".epi",core_file)[i],
         read.table(paste0(epi_path,epi_files[i]),
                    col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele")
         )    
  )
  assign(gsub(".class.txt",".filtered",core_file)[i],
         merge(get(CDS_name[i]),
               get(gsub(".class.txt",".epi",core_file)[i]),
               by = "gene", all.x = T)
  )
  assign(gsub(".class.txt",".aggre",core_file)[i],
         aggregate(epiallele~pan,
                   get(gsub(".class.txt",".filtered",core_file)[i]), 
                   paste %>% unlist(use.names = F)))
} 
#Make the matrix with epiallele information
core_pan_ID <- read.table("/Users/x/Desktop/Data/pan_gene/core_panID.txt",
                          col.names = "pan")
epi_matrix <- merge(core_pan_ID,select(B73.aggre,c("pan","epiallele")), by = "pan", all.x = T)
for (i in 2:26) {
  epi_matrix <- merge(epi_matrix,
                      get(gsub(".class.txt",".aggre",core_file)[i]),
                      by = "pan", all.x = T)
}
names(epi_matrix) <- c("pan",gsub(".class.txt","",core_file))

#Define the epiallele status for pangenes
epi_matrix$UM = rowSums(epi_matrix[,2:27]=="UM")
epi_matrix$gbM = rowSums(epi_matrix[,2:27]=="gbM")
epi_matrix$teM = rowSums(epi_matrix[,2:27]=="teM")

epi_matrix$status ="unknown"
epi_matrix$status[epi_matrix$UM>=2&epi_matrix$gbM==0&epi_matrix$teM==0] <- "UM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM>=2&epi_matrix$teM==0] <- "gbM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM==0&epi_matrix$teM>=2] <- "teM"
epi_matrix$status[epi_matrix$UM>=2&epi_matrix$gbM>=2&epi_matrix$teM==0] <- "UM_gbM"
epi_matrix$status[epi_matrix$UM>=2&epi_matrix$gbM==0&epi_matrix$teM>=2] <- "UM_teM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM>=2&epi_matrix$teM>=2] <- "gbM_teM"
epi_matrix$status[epi_matrix$UM>=2&epi_matrix$gbM>=2&epi_matrix$teM>=2] <- "UM_gbM_teM"

#import the copy number status for
pan_copy_1_1 <-  cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/1-1_pangene.txt",
           header = T), "1-1")
pan_copy_1_N <- cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/1-N_pangene.txt",
                 header = T), "1-N")
pan_copy_N_N <- cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/N-N_pangene.txt",
                 header = T), "N-N")
names(pan_copy_1_1) <-names(pan_copy_1_N) <- names(pan_copy_N_N)<-  c("pan","copy")

##Calculate the number of pangenes in each category
pan_copy <- rbind(pan_copy_1_1,pan_copy_1_N, pan_copy_N_N)
epi_matrix <- merge(epi_matrix,pan_copy, by = "pan", all.x = T)
table(epi_matrix$status,epi_matrix$copy)
sum(table(epi_matrix$status,epi_matrix$copy)[c(1:7),2])
sum(table(epi_matrix$status,epi_matrix$copy)[c(1,3,4),2])
epi_matrix_final <- epi_matrix %>%filter(copy == "1-N") %>% filter(status %in% c("UM","gbM","teM"))
table(epi_matrix_final$status)

#Write out the epiallele matrix for 1-N pangenes
write_rds(epi_matrix_final,"/Users/x/Desktop/epi_matrix_final.rda")

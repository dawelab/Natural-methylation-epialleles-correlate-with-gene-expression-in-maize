library(tidyverse)
library(matrixStats)
#perform the filtering by cumulative CDS lengths
pan_matrix <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv")
(NAM_CDS <- names(pan_matrix)[4:29])
github_path <- "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Methylation_Expression/CDS/"
NAM_CDS_list <- c(paste0("Zm-",NAM[1],"-REFERENCE-NAM-5.0.1.canon.cds_length.bed"),
           paste0("Zm-",NAM[2:26],"-REFERENCE-NAM-1.0.1.canon.cds_length.bed"))
#Import CDS length for each CDS element
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".CDS"),
         read.table(paste0(github_path,NAM_CDS_list[i]),
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
coregene_path = "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/"
core_file <- paste0(NAM_CDS,".class.txt")
for (i in 1:length(core_file)) {
  assign(gsub(".txt","",core_file[i]),read.table(paste0(coregene_path,core_file[i]),
                                                 sep = ",", header = T) %>%
           filter(class == "Core Gene")
  )
}
#cbind(core_file,NAM_CDS_list)
for (i in 1:length(NAM_CDS)) {
  assign(paste0(NAM_CDS[i],".preCDS"),
         merge(get(gsub(".txt","",core_file[i])),
               get(paste0(NAM_CDS[i],".comCDS")),
               by="gene", all.x = T) %>% select(c("pan","gene","length")) 
         
  )
}
#Read in the CDS matrix to extract the pangene and its median
CDS_matrix <- read.table("/Users/x/Desktop/NAM_CDS_matrix.txt",
                         header = T)[,c("pan","median")]
#Combine the median of CDS length together for filtering genes with 10% difference in CDS length
CDS_name <- gsub(".class.txt",".panCDS",core_file)
for (i in 1:26){
  assign(CDS_name[i],
         merge(CDS_matrix,get(paste0(NAM_CDS[i],".preCDS")),by="pan", all.x = T) %>% 
           filter(abs(length-median)/median <= 0.1)) %>% filter(!is.na(length))
}
#Subsetting the epialleles for filtered gene IDs
#Read in the methylation files
epi_path <- "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/epiallele/"
epi_files  <- c(paste0("Zm-",NAM[1],"-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt"), #Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt
                  paste0("Zm-",NAM[2:26],"-REFERENCE-NAM-1.0.1.canon.gene.gene.mtr.ID.type.txt"))
#
#cbind(epi_files,core_file)
for (i in 1:26) {
  assign(gsub(".class.txt",".epi",core_file)[i],
         read.table(paste0(epi_path,epi_files[i]),
                    col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele")
         )    
  )
}
for (i in 1:26) {
  assign(gsub(".class.txt",".filtered",core_file)[i],
         merge(get(CDS_name[i]),
               get(gsub(".class.txt",".epi",core_file)[i]),
               by = "gene", all.x = T)
  )
}
for (i in 1:26) {
  assign(gsub(".class.txt",".aggre",core_file)[i],
         aggregate(epiallele~pan,
                   get(gsub(".class.txt",".filtered",core_file)[i]), 
                   paste %>% unlist(use.names = F)))
} 
#Make the matrix with epiallele information
core_pan_ID <- read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/core_panID.txt",
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
epi_matrix$N =  rowSums(epi_matrix[,c("UM","gbM","teM")])
epi_matrix$status ="unknown"
epi_matrix$status[epi_matrix$UM!=0&epi_matrix$gbM==0&epi_matrix$teM==0] <- "UM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM!=0&epi_matrix$teM==0] <- "gbM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM==0&epi_matrix$teM!=0] <- "teM"
epi_matrix$status[epi_matrix$UM!=0&epi_matrix$gbM!=0&epi_matrix$teM==0] <- "UM_gbM"
epi_matrix$status[epi_matrix$UM!=0&epi_matrix$gbM==0&epi_matrix$teM!=0] <- "UM_teM"
epi_matrix$status[epi_matrix$UM==0&epi_matrix$gbM!=0&epi_matrix$teM!=0] <- "gbM_teM"
epi_matrix$status[epi_matrix$UM!=0&epi_matrix$gbM!=0&epi_matrix$teM!=0] <- "UM_gbM_teM"

#import the copy number status for
pan_copy_1_1 <-  cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/1-1_pangene.txt",
                                  header = T), "1-1")
pan_copy_1_N <- cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/1-N_pangene.txt",
                                 header = T), "1-N")
pan_copy_N_N <- cbind(read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Epiallele%20stability/N-N_pangene.txt",
                                 header = T), "N-N")
names(pan_copy_1_1) <-names(pan_copy_1_N) <- names(pan_copy_N_N)<-  c("pan","copy","duplicate")

##Calculate the number of pangenes in each category
pan_copy <- rbind(pan_copy_1_1,pan_copy_1_N, pan_copy_N_N)
epi_matrix <- merge(epi_matrix,pan_copy, by = "pan", all.x = T)

#For pangenes with 2 exist as singleton and as least one exist as duplicate
epi_matrix %>% filter(duplicate == "1-N") %>% filter(N >=2) %>% dim
#singletons were present in more than one epiallele state in 1-N pangenes
epi_matrix %>% filter(duplicate == "1-N") %>% filter(N >=2) %>% filter(!(status %in% c("UM","gbM","teM"))) %>% dim
#represented by singletons of one epiallele type in at least two genomes
epi_matrix %>% filter(duplicate == "1-N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>% dim
#The distribution of the maximum copy number
epi_matrix %>% filter(duplicate == "1-N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(copy) %>% table()

#The distribution of epiallele in 1_N distribution
epi_matrix %>% filter(duplicate == "1-N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(status) %>% table()

#The joint distribution of epiallele in 1_N % 1_1 distribution as stable panegen
(joint11_1N <- epi_matrix  %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(status,duplicate) %>% table())
joint11_1N[,2]/joint11_1N[,1]
0.5926296/0.2650200

#The joint distribution of epiallele in 1_N % 1_1 distribution as all pangene
 epi_matrix  %>% filter(N >=2) %>%  select(status,duplicate) %>% table()



#Select the final 1-N matrix and export as rds format
epi_matrix_final <- epi_matrix %>%filter(duplicate == "1-N") %>% filter(status %in% c("UM","gbM","teM")) %>% filter(N >=2)

#Write out the epiallele matrix for 1-N pangenes
write_rds(epi_matrix_final,"/Users/x/Desktop/epi_matrix_final.rda")

#Select the final 1-1 stable Matrice and export as rds format
epi_matrix_1_1_UM_gbM <- epi_matrix %>%filter(duplicate == "1-1" & UM >=2 & gbM >=2) %>% filter(status %in% c("UM_gbM")) 
epi_matrix_1_1_UM_teM <- epi_matrix %>%filter(duplicate == "1-1" & UM >=2 & teM >=2) %>% filter(status %in% c("UM_teM")) 
epi_matrix_1_1_gbM_teM <- epi_matrix %>%filter(duplicate == "1-1" & gbM >=2 & teM >=2) %>% filter(status %in% c("gbM_teM")) 

write_rds(epi_matrix_1_1_UM_gbM,"/Users/x/Desktop/epi_matrix_1_1_gbM_teM.rda")
write_rds(epi_matrix_1_1_UM_teM,"/Users/x/Desktop/epi_matrix_1_1_UM_teM.rda")
write_rds(epi_matrix_1_1_gbM_teM,"/Users/x/Desktop/epi_matrix_1_1_gbM_teM.rda")

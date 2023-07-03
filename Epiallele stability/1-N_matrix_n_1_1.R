library(tidyverse)
library(matrixStats)
#perform the filtering by cumulative CDS lengths
path <- "/Users/x/Desktop/Data/combine/"
files <- list.files(path)
for (i in 1:26) {
  assign(
    gsub(".csv","",files[i]),
    read.csv(paste0(path,files[i]))
  )
}
NAM <- gsub(".all.csv","",files)
pan_coreID <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv") %>% 
  filter(class == "Core Gene") %>% select("Pan_gene_ID")
names(pan_coreID) <- "pan"

#Read in the CDS matrix to extract the pangene and its median
CDS_matrix <- read.csv("/Users/x/Desktop/CDS_median_matrix.csv",
                         header = T)[,c("pan","median")]
#Combine the median of CDS length together for filtering genes with 10% difference in CDS length
CDS_name <- gsub(".all.csv",".panCDS",files)
for (i in 1:26){
  assign(CDS_name[i],
         merge(CDS_matrix,get(gsub(".csv","",files[i])),by="pan", all.x = T) %>% 
           filter(abs(length-median)/median <= 0.1)) %>% filter(!is.na(length))
}

###Select pan and epiallele info to generate the epiallele matrix and aggregate
epi_matrix <- merge(pan_coreID,
                    aggregate(epiallele~pan,
                              get(CDS_name[1]) %>% select(c("pan","epiallele")),
                              FUN = paste0%>%unlist(use.names = F)),
                    by = "pan", all.x = T)

for (i in 2:26) {
  epi_matrix<-  merge(epi_matrix,
        aggregate(epiallele~pan,
                  get(CDS_name[i]) %>% select(c("pan","epiallele")),
                  FUN = paste0%>%unlist(use.names = F)),
        by = "pan", all.x = T)
}

names(epi_matrix) <- c("pan",NAM)

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


pan_copy <- copy_matrix %>% select(c("pan","status","max"))
names(pan_copy) <- c("pan","duplicate","copy")
epi_matrix <- merge(epi_matrix,pan_copy, by = "pan", all.x = T)

#For pangenes with 2 exist as singleton and as least one exist as duplicate
epi_matrix %>% filter(duplicate == "1_N") %>% filter(N >=2) %>% dim
#singletons were present in more than one epiallele state in 1-N pangenes
epi_matrix %>% filter(duplicate == "1_N") %>% filter(N >=2) %>% filter(!(status %in% c("UM","gbM","teM"))) %>% dim
#represented by singletons of one epiallele type in at least two genomes
epi_matrix %>% filter(duplicate == "1_N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>% dim
#The distribution of the maximum copy number
epi_matrix %>% filter(duplicate == "1_N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(copy) %>% table()

#The distribution of epiallele in 1_N distribution
epi_matrix %>% filter(duplicate == "1_N") %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(status) %>% table()

#The joint distribution of epiallele in 1_N % 1_1 distribution as stable panegen
(joint11_1N <- epi_matrix  %>% filter(N >=2) %>% filter((status %in% c("UM","gbM","teM"))) %>%  select(status,duplicate) %>% table())
joint11_1N[,2]/joint11_1N[,1]
0.5926296/0.2651869

#The joint distribution of epiallele in 1_N % 1_1 distribution as all pangene
epi_matrix  %>% filter(N >=2) %>%  select(status,duplicate) %>% table()



#Select the final 1-N matrix and export as rds format
epi_matrix_final <- epi_matrix %>%filter(duplicate == "1-N") %>% filter(status %in% c("UM","gbM","teM")) %>% filter(N >=2)

#Write out the epiallele matrix for 1-N pangenes
write_rds(epi_matrix_final,"/Users/x/Desktop/epi_matrix_final.rds")

#Select the final 1-1 stable Matrice and export as rds format
epi_matrix_1_1_UM_gbM <- epi_matrix %>%filter(duplicate == "1_1" & UM >=2 & gbM >=2) %>% filter(status %in% c("UM_gbM")) 
epi_matrix_1_1_UM_teM <- epi_matrix %>%filter(duplicate == "1_1" & UM >=2 & teM >=2) %>% filter(status %in% c("UM_teM")) 
epi_matrix_1_1_gbM_teM <- epi_matrix %>%filter(duplicate == "1_1" & gbM >=2 & teM >=2) %>% filter(status %in% c("gbM_teM")) 
epi_matrix_1_1_gbM_teM <- epi_matrix %>%filter(duplicate == "1_1" & UM >=2 & gbM >=2 & teM >=2) %>% filter(status %in% c("UM_gbM_teM")) 

write_rds(epi_matrix_1_1_UM_gbM,"/Users/x/Desktop/epi_matrix_1_1_gbM_teM.rds")
write_rds(epi_matrix_1_1_UM_teM,"/Users/x/Desktop/epi_matrix_1_1_UM_teM.rds")
write_rds(epi_matrix_1_1_gbM_teM,"/Users/x/Desktop/epi_matrix_1_1_gbM_teM.rds")

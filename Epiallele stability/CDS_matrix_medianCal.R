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

#Generate the CDS matrix for median calculation
for (i in 1:26) {
  assign(
    paste0(NAM[i],".CDS"),
    get(gsub(".csv","",files[i])) %>% select(c("pan","length")) %>% unique()
  )
}
##Aggregate the multi-copy gene into one cell
for (i in 1:26) {
  assign(
    paste0(NAM[i],".CDS.agg"),
    aggregate(length~pan,get(paste0(NAM[i],".CDS")),FUN = paste0 %>% unlist)
  )
}
##
CDS_matrix <- merge(pan_coreID,B73.CDS.agg,by="pan",all.x = T)
for (i in 2:26) {
  CDS_matrix <- merge(CDS_matrix,
                      get(paste0(NAM[i],".CDS.agg")),
                      by = "pan",
                      all.x = T)
}
names(CDS_matrix)  <- c("pan",NAM)
##Make the multicopy gene as NA
sum(CDS_matrix$pan!=copy_matrix$pan)
sum(names(CDS_matrix)[1:27]!=names(copy_matrix)[1:27])
CDS_matrix_ij <- CDS_matrix[,2:27]
copy_matrix_ij <- copy_matrix[,2:27]
CDS_matrix_ij[copy_matrix_ij>1] <- NA
CDS_matrix_ij_num <-  as.matrix(CDS_matrix_ij) %>% as.numeric() %>% matrix(ncol=26)

pan_median <- rowMedians(CDS_matrix_ij_num,na.rm = T)
CDS_matrix_ij_num = as.data.frame(CDS_matrix_ij_num)
CDS_matrix_ij_num$pan = CDS_matrix$pan
CDS_matrix_ij_num$median = pan_median
names(CDS_matrix_ij_num) <- c(names(CDS_matrix)[2:27],"pan","median")
CDS_matrix_ij_num <- CDS_matrix_ij_num[,c(names(CDS_matrix),"median")]

write.csv(
  CDS_matrix_ij_num,"/Users/x/Desktop/CDS_median_matrix.csv",
  row.names = F
)

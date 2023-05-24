setwd("/Users/x/Desktop/Data/matrix/CDS/")
NAM_CDS_list = list.files()
NAM_CDS <- gsub("-REFERENCE-NAM-1.0.1.canon.cds_length.bed","",NAM_CDS_list)
NAM_CDS <- gsub("-REFERENCE-NAM-5.0.1.canon.cds_length.bed","",NAM_CDS)
NAM_CDS <- gsub("Zm-","",NAM_CDS)
NAM_CDS[26] <- "Tzi8"
#Import CDS length for each CDS elementg
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
coregene_path = "/Users/x/Desktop/Data/core/"
core_file <- list.files(coregene_path)
for (i in 1:length(core_file)) {
  assign(core_file[i],read_table(paste0(coregene_path,core_file[i]), col_names = 
                                   c("chr","start","end","strand","gene","pan","copy","duplicate"))[,5:8])
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
names(NAM_CDS_matrix) <- c("pan",gsub(".core","",core_file))
#Generate a copy matrix for NAM founders
copy_list <- gsub("core","copy",core_file)
for (i in 1:length(copy_list)) {
  assign(copy_list[i],
         unique(select(get(core_file[i]),c("pan","copy")))
         )
}
NAM_copy_matrix <- merge(get(copy_list[1]),
                         get(copy_list[2]),
                         by="pan",
                         all=T) 
for (i in 3:26) {
  NAM_copy_matrix <- merge(NAM_copy_matrix,
                           get(copy_list[2]),
                           by="pan",
                           all=T)
}
names(NAM_copy_matrix) <- c("pan",gsub(".core","",core_file))
#Compare index mCG matrix and copy number matrix
sum(names(NAM_copy_matrix) != names(NAM_CDS_matrix))
sum(NAM_copy_matrix$pan!=NAM_CDS_matrix$pan)
#Transform into double type
NAM_CDS_matrix[NAM_CDS_matrix=="NULL"] <- NA
for (i in 1:26) {
      for (j in 1:dim(NAM_CDS_matrix)) {
        if("," %in% NAM_CDS_matrix[j,i+1]) NAM_CDS_matrix[j,i] = NA
      }  
}

NAM_CDS_matrix <- as.character( NAM_CDS_matrix)
 NAM_CDS_matrix["," %in% NAM_CDS_matrix]


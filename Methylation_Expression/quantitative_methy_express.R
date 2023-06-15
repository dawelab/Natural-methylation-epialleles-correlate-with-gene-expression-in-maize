library(tidyverse)
#Match the pangene and gene name to create the pan matrix for the CG/CHG methylation levels
#1.Read all core  gene files
coregene_path = "/Users/x/Desktop/Data/core/"
core_file <- list.files(coregene_path)
for (i in 1:length(core_file)) {
  assign(core_file[i],read_table(paste0(coregene_path,core_file[i]), col_names = 
             c("chr","start","end","strand","gene","pan","copy","duplicate"))[,5:8])
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

#mask the ambiguous and teM genes' mCG value
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

#get rid of the epiallele colum
for (i in 1:length(mCG_list)) {
  assign(mCG_list[i],
       select(get(pan_methy_names[i]),
              c("pan","mCG")
              )
        )

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
NAM_CDS_matrix_int <- read.table("/Users/x/Desktop/NAM_CDS_matrix.txt")
mCG_pan_list = data.frame(pan=pan_mCG$pan)
CDS_matrix <- merge(mCG_pan_list,
                    NAM_CDS_matrix_int,
                    by = "pan",
                    all.x = T)
CDS_matrix<- CDS_matrix[,-28]
sum(names(pan_mCG) != names(CDS_matrix))
sum(pan_mCG$pan != CDS_matrix$pan)
pan_mCG[is.na(CDS_matrix)] <- NA
# Make expression data
express <- "/Users/x/Desktop/Data/expression/TPM/"
express_file <- list.files(express)
express_list <- gsub(".txt","",express_file)
for (i in 1:26) {
  assign(express_list[i],
         select(read.table(paste0(express,express_file[i]),
                    header = T),
                c("gene","middle"))
         )
}
#match pangene & gene ID
for (i in 1:26) {
  assign(paste0(express_list[i],".middle"),
       select(merge(get(express_list[i]),
             get(core_file[i]),
             by = "gene",
             all.y = T),
             c("pan","middle"))
)
}
#aggregate by pangene
for (i in 1:26) {
  assign(paste0(express_list[i],".pre"),
       aggregate(middle~pan,
                 data = get(paste0(express_list[i],".middle")),
                 FUN = paste0%>%unlist
                 )
)
  ##aggregate by pangene  for copy number
for (i in 1:26) {
  assign(paste0(express_list[i],".copy"),
       aggregate(middle~pan,
                 data = get(paste0(express_list[i],".middle")),
                 FUN = length
                 )
)
}
#Create the middle expression matrix
middle_matrix <- merge(get(paste0(express_list[1],".pre")),
                       get(paste0(express_list[2],".pre")),
                       by="pan",
                       all = T)
for (i in 3:26) {
  middle_matrix <- merge(middle_matrix,
                         get(paste0(express_list[i],".pre")),
                         by="pan",
                         all = T
                         )
}
names(middle_matrix) <- c("pan",gsub(".core","",core_file))

#Create a matrix for matrix copy number
#Create the middle expression matrix
middle_copy <- merge(get(paste0(express_list[1],".copy")),
                       get(paste0(express_list[2],".copy")),
                       by="pan",
                       all = T)
for (i in 3:26) {
  middle_copy <- merge(middle_copy,
                         get(paste0(express_list[i],".copy")),
                         by="pan",
                         all = T
                         )
}
names(middle_copy) <- c("pan",gsub(".core","",core_file))

#Mask the genes wtih big CDS difference
middle <- merge(mCG_pan_list,
      middle_matrix,
      by="pan",
      all.x = T)

middle[is.na(pan_mCG)] <- NA
middle <- as.matrix(middle)
middle <- as.numeric(middle)
middle <- matrix(middle, ncol)
cor_list <- c()
for (i in c(1:29,31:dim(middle)[1])){
  x = as.numeric(as.matrix(unlist(middle[i,2:27])))
  y = as.numeric(as.matrix(unlist(pan_mCG[i,2:27])))
  if(sum(!is.na(x)) >=3) {
  stat = cor.test(x,y)$estimate
  cor_list <- c(cor_list,stat)
  }
  if(sum(!is.na(x)) < 3) {
  stat = NA
  cor_list <- c(cor_list,stat)
  }
}

library(tidyverse)
library(matrixStats)
#Match the pangene and gene name to create the pan matrix for the CG/CHG methylation levels
#1.Read all core  gene files
pan_matrix <- read.csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv")
(NAM <- names(pan_matrix)[4:29])

coregene_path = "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/"
  core_file <- paste0(NAM,".class.txt")

for (i in 1:length(core_file)) {
  assign(core_file[i],read.table(paste0(coregene_path,core_file[i]), 
                                 header = T, sep = "," ) %>% filter(class == "Core Gene" & copy == 1))
}
#2. Read all methylation level files
methy_path = "https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/epiallele/"
methy_file = c("Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt",
                 paste0("Zm-",NAM[2:26],"-REFERENCE-NAM-1.0.1.canon.gene.gene.mtr.ID.type.txt"))
  
methy_names = gsub(".class.txt","methy",core_file)
for (i in 1:length(methy_file)) {
    assign(methy_names[i],
           read.table(paste0(methy_path,methy_file[i]),
                      col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele")) %>% 
             filter(cCG >= 40 & cCHG >= 40)) 
  }
#3. Use merge to put together pangene names and gene names
mCG_names = gsub(".class.txt","mCG",core_file)
for (i in 1:length(pan_methy_names)) {
    assign(mCG_names[i],
           merge(get(core_file[i]),
                 get(methy_names[i]),
                 by = "gene") %>% filter(mCHG <= 0.05) %>% select(c("mCG","pan")) )
}


mCHG_names = gsub(".class.txt","mCHG",core_file)
for (i in 1:length(pan_methy_names)) {
  assign(mCHG_names[i],
         merge(get(core_file[i]),
               get(methy_names[i]),
               by = "gene") %>% select(c("mCHG","pan")) )
}


#Read in the core gene panID
core_panID <- read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/core_panID.txt",
                         col.names = "pan")
#Make the mCG matrix
pan_mCG <- merge(core_panID,
                  get(mCG_names[1]),
                  by = "pan",
                  all.x = T)
for (i in 2:length(mCG_names)) {
  pan_mCG <- merge(pan_mCG,
                    get(mCG_names[i]),
                    by="pan",
                    all.x=T)
}
names(pan_mCG) <- c("pan",gsub(".class.txt","",core_file))
#Make the mCHG matrix
pan_mCHG <- merge(core_panID,
                     get(mCHG_names[1]),
                     by = "pan",
                     all.x = T)
for (i in 2:length(mCHG_names)) {
      pan_mCHG <- merge(pan_mCHG,
                       get(mCHG_names[i]),
                       by="pan",
                       all.x=T)
  }
names(pan_mCHG) <- c("pan",gsub(".class.txt","",core_file))
#Filter out the mCG values with CDS length in big difference
CDS_matrix <- read.table("/Users/x/Desktop/NAM_CDS_matrix.txt", header = T)

CDS_matrix<- CDS_matrix[,-28]
sum(names(pan_mCHG) != names(CDS_matrix))
sum(pan_mCHG$pan != CDS_matrix$pan)
pan_mCHG[is.na(CDS_matrix)] <- NA
pan_mCG[is.na(CDS_matrix)] <- NA

pan_mCG_diff <- rowMaxs(as.matrix(pan_mCG[,2:27]),na.rm = T) - rowMins(as.matrix(pan_mCG[,2:27]),na.rm = T)

pan_mCG_mat <- pan_mCG[pan_mCG_diff>0.2,] 
pan_mCHG_mat <- pan_mCHG[pan_mCHG_diff>0.2,] 

write.csv(pan_mCG_mat,"/Users/x/Desktop/Data/matrix/pan_mCG_mat.csv",quote = F)
write.csv(pan_mCHG_mat,"/Users/x/Desktop/Data/matrix/pan_mCHG_mat.csv",quote = F)

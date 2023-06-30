##Methylation files
class_path = "/Users/x/Desktop/Data/core_gene/"
class_file <- list.files(class_path)[-c(11,12)]

for (i in 1:length(class_file)) {
  assign(class_file[i],read.table(paste0(class_path,class_file[i]), 
                                 header = T, sep = "," ) )
}
methy_path = "/Users/x/Desktop/Data/methylation/cgchgmtr/"
methy_file = list.files(methy_path)
methy_names = gsub("class.txt","methy",class_file)

for (i in 1:length(methy_file)) {
  assign(methy_names[i],
         read.table(paste0(methy_path,methy_file[i]),
                    col.names = c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele")))
}
# Use merge to put together pangene names and gene names
#cbind(methy_file,class_file)


merge_names = gsub("class.txt","merge",class_file)
for (i in 1:length(mCG_names)) {
  assign(merge_names[i],
         merge(get(class_file[i]),
               get(methy_names[i]),
               by = "gene", all.y = T)[,c("chr","start","end","strand","mCG","cCG","mCHG","cCHG","gene","epiallele","pan","class","copy")])
}

#Merge with tpm files
cbind(merge_names,combine_names,tpm_file)
tpm_path <- "/Users/x/Desktop/Data/expression/TPM/"
tpm_file <- list.files(tpm_path)
for (i in 1:26) {
  assign(tpm_file[i],
         read.table(paste0(tpm_path,tpm_file[i]),header = T)[,-c(1:3)])
}

combine_names = gsub("class.txt","combine",class_file)
for (i in 1:26) {
  assign(combine_names[i],
         merge( get(merge_names[i]),
           get(tpm_file[i]),
           by="gene", all.x = T
         ))
}

##Combine with the cumulative CDS length
CDS_path  <- "/Users/x/Desktop/Data/matrix/CDS/"
CDS.file <- list.files(CDS_path)[-1]
CDS.name <- gsub("class.txt","CDS",class_file)

for (i in 1:26) {
  assign(CDS.name[i],
         read.table(paste0(CDS_path,CDS.file[i]),
                    col.names = c("gene", "length")))
  assign(paste0(CDS.name[i],".com"),
         aggregate(length~gene,get(CDS.name[i]),sum)
         )
}

#Combine all the information into one files
all_names = gsub("class.txt","all",class_file)
for (i in 1:26) {
  assign(all_names[i],
         merge( get(combine_names[i]),
                get(paste0(CDS.name[i],".com")),
                by="gene", all.x = T
         ))
}
out <- "/Users/x/Desktop/Data/combine/"
for (i in 1:26) {
  write.csv(get(all_names[i]),
            paste0(out,all_names[i],".csv"), quote = F, row.names = F)
}

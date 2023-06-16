setwd("/Users/x/Desktop/Data/core_gene")
melt = function(df,NAM){
  gene = NULL
  for (i in 1:(dim(df)[2]-2)){
    temp = as.character(df[,i])
    gene = c(gene,temp)
  }
  
  pan = rep(df$pan,dim(df)[2]-2)
  class = rep(df$pangeneclass,dim(df)[2]-2)
  
  df2 = data.frame(gene,pan,class)
  df3 = df2[complete.cases(df2),]
  return(write.csv(df3,paste0(NAM,".class"), sep = "\t", quote = F, row.names = F))
}
samples = list.files()
sample_name <- gsub(".class.txt","",samples)
for (i in 1:26){
  data = read.csv(samples[i],sep = ",", header = T)
  melt(data,sample_name[i])
}

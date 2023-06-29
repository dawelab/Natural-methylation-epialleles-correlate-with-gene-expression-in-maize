setwd("/Users/x/Desktop/Data/core_gene")
melt = function(df,NAM){
  
  #melt all the gene into list and remove NA and organize into dataframe
  gene = NULL
  for (i in 1:(dim(df)[2]-2)){
    temp = as.character(df[,i])
    gene = c(gene,temp)
  }
  
  pan = rep(df$pan,dim(df)[2]-2)
  class = rep(df$pangeneclass,dim(df)[2]-2)

  df2 = data.frame(gene,pan,class)
  df3 = df2[complete.cases(df2),]
  copy = aggregate(gene~pan,data=df3, FUN= length)
  names(copy)[2] <- "copy"
  df4 =  merge(df3,copy,by="pan")
  return(write.csv(df4,paste0(NAM,".class"), sep = ",", quote = F, row.names = F))
}
samples = list.files()[-c(11:12)]
sample_name <- gsub(".class.txt","",samples)
for (i in 1:26){
  data =  read.csv(samples[i],sep = ",", header = T)
  melt(data,sample_name[i])
}

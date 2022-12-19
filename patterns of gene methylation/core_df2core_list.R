setwd("/Users/x/Desktop/Data/pan_gene/version3")
melt = function(df,NAM){
  gene = NULL
  for (i in 1:(dim(df)[2]-3)){
    temp = as.character(df[,i])
    gene = c(gene,temp)
  }
  
  copy = rep(df$copy,dim(df)[2]-3)
  pan = rep(df$pan,dim(df)[2]-3)
  class = rep(df$class,dim(df)[2]-3)
  
  df2 = data.frame(gene,pan,copy,class)
  df2$duplicate = "Y"
  df2$duplicate[df2$copy==1] <- "N"
  df3 = df2[complete.cases(df2),]
  
  return(write.csv(df3,NAM, sep = "\t", quote = F, row.names = F))
}
samples = list.files("/Users/x/Desktop/Data/pan_gene/version3")[-24]
for (i in 1:26){
  data = read.table("B97.core",sep = ",", header = T)
  melt(data,"B97")
}

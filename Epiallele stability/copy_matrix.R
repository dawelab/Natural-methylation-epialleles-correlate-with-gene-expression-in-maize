library(matrixStats)
#Read in the files with pangene class information 
coregene_path = "/Users/x/Desktop/Data/core_gene/"

core_file <- list.files(coregene_path)
for (i in 1:length(core_file)) {
  assign(gsub(".txt","",core_file[i]),
         read.table(paste0(coregene_path,core_file[i]), 
                    sep = ",", header = T ) %>% filter(class == "Core Gene"))
}
#Generate the matrix with copy number information
copy_name <- gsub(".class.txt",".copy", core_file)
for (i in 1:length(core_file)) {
  assign(copy_name[i],
         aggregate(gene~pan, 
                   get(core_file[i]),
                   length)
  )  
}
#Loop the copy information into a matrix
#Read in the core gene panID
core_pan <- read.table("/Users/x/Desktop/Data/pan_gene/core_panID.txt",
                       col.names = "pan")
copy_matrix <- merge(core_pan,
                     B73.copy,
                     by = "pan",
                     all.x = T)
for (i in 2:26) {
  copy_matrix <- merge(copy_matrix,
                       get(copy_name[i]),
                       by = "pan",
                       all.x = T)
}
names(copy_matrix) <- c("pan",gsub(".copy","",copy_name))
#Calculate the minimum and maximum copy number for each pangene
copy_matrix$min <- rowMins(as.matrix(copy_matrix[,2:27]),na.rm = T)
copy_matrix$max <- rowMaxs(as.matrix(copy_matrix[,2:27]),na.rm = T)
#1-1 singleton pangenes
copy_matrix$class <- "unknown"
copy_matrix$class[copy_matrix$min==1 & copy_matrix$max==1] <- "1-1"
copy_matrix$class[copy_matrix$min==1 & copy_matrix$max>1] <- "1-N"
copy_matrix$class[copy_matrix$min>1 & copy_matrix$max>1] <- "N-N"
table(copy_matrix$class) %>% sum
write.table(copy_matrix[copy_matrix$class == "1-1", c("pan","max")], 
            "/Users/x/Desktop/Data/core_gene/1-1_pangene.txt",
            quote = F,
            row.names = F, col.names = c("pan","max"))
write.table(copy_matrix[copy_matrix$class == "1-N", c("pan","max")], 
            "/Users/x/Desktop/Data/core_gene/1-N_pangene.txt",
            quote = F,
            row.names = F, col.names = c("pan","max"))
write.table(copy_matrix[copy_matrix$class == "N-N", c("pan","max")], 
            "/Users/x/Desktop/Data/core_gene/N-N_pangene.txt",
            quote = F,
            row.names = F, col.names = c("pan","max"))

#Read in the  table with the gene list syntenic to sorhgam
sorhgum <- read.table("/Users/x/Desktop/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_Sb_subgenomes_complete_longest.bed",
                      sep = "\t")[,4]%>% unique()  %>% as.data.frame() 
names(sorhgum) <- "gene"
sorhgum$synt <- "Y"
#To examine if there is any na in the provided sorhgam table 
any(is.na(sorhgum))

#Read in the table with B73 pangene class & epiallele information
#Here is the code for calculating B73 core gene set
B73 <- read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/epiallele/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt")[,9:10]
names(B73) <- c("gene", "epiallele")  
B73_class <- read.table("https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/patterns%20of%20gene%20methylation/pangene_class/B73.class.txt",
                        sep = ",", header = T)
B73_core <- merge(B73,B73_class, by = "gene") %>% filter(class == "Core Gene")
B73_pan <- merge(B73,B73_class, by = "gene") 

#Here is the code for calculating B73 all gene set
combi_all <- merge(B73_pan,sorhgum, by = "gene", all.x = T)
combi_all$synt[is.na(combi_all$synt)] <- "N"
(count_all <- table(combi_all$epiallele,combi_all$synt)[c(4,2,3),2:1])
rowSums(count_all)
count_all/rowSums(count_all)

#Here is the code for calculating B73 core gene set
combi_core <- merge(B73_core,sorhgum, by = "gene", all.x = T)
combi_core$synt[is.na(combi_core$synt)] <- "N"
(count_core <- table(combi_core$epiallele,combi_core$synt)[c(4,2,3),2:1])
count_core/rowSums(count_core)

##ploting
df <- data.frame(proportion = c(count_all/rowSums(count_all),count_core/rowSums(count_core)),
                 epiallele = rep(c("UM","gbM","teM"),4),
                 Syntenic = rep(rep(c("Y","N"),each=3),2),
                 gene = rep(c("all","core"),each = 6)
                 )
df$epiallele = factor(df$epiallele, levels = c("UM","gbM","teM"))
df %>% filter(Syntenic=="Y") %>% ggplot(aes(x=epiallele,y=proportion, fill = gene)) +
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           alpha = 0.75,
           width = 0.35) +
  theme_bw() +
  ylab("proportion with sorghum synteny") + xlab("") +
  scale_fill_manual(values=c("#929591","#000000"))
#Replace the all with grey color and core with black color
#legend title - gene type 



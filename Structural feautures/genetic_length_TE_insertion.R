library(tidyverse)
library(extrafont)
#Read in different element length
CDS=read.table("/Users/x/Desktop/Data/TE/Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS.length.bed")
exon=read.table("/Users/x/Desktop/Data/TE/Zm-B73-REFERENCE-NAM-5.0.1.canon.exon.length.bed")
UTR = read.table("/Users/x/Desktop/Data/TE/Zm-B73-REFERENCE-NAM-5.0.1.canon.UTR.length.bed")
gene = read.table("/Users/x/Desktop/Data/TE/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.length")
TE = cbind(read.table("/Users/x/Desktop/Data/TE/B73_mergingTE.interset.bed")[,c(5,9)],"TE")

#Read in core list in B73
core_list2  = data.frame(gene=read.table("/Users/x/Desktop/Data/core/B73.core")[,5])

names(exon)  <- names(UTR) <- names(CDS) <- names(TE) <- c("gene","length","element")
names(gene) <- c("gene","length") 

B73.core = read.table("/Users/x/Desktop/Data/pan_gene/version3/B73.class.txt.txt",
                      sep=",",header = T)
B73_core_list = data.frame(gene=B73.core$gene[B73.core$class=="Core Gene"])
for (i in 1:dim(B73_core_list)[1]) {
  B73_core_list[i,1] = gsub("^gmap_ID=chr",NA,B73_core_list[i,1])
  
}
B73_core_list = B73_core_list$gene[!is.na(B73_core_list$gene)]
B73_core_list = data.frame(gene=B73_core_list)
B73_epiallele = read.table("/Users/x/Desktop/Data/methylation/cgchgmtr/Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt")[,c(9,10)]
names(B73_epiallele) = c("gene","epiallele")
B73_core_epiallele = merge(B73_core_list,B73_epiallele,by="gene",all.x = T)
B73_core_epiallele = B73_core_epiallele[B73_core_epiallele$epiallele!="ambiguous",]
UTR_combine = cbind(aggregate(UTR$length,by = list(UTR$gene),sum),"UTR")
exon_combine = cbind(aggregate(exon$length,by = list(exon$gene),sum),"exon")
CDS_combine = cbind(aggregate(CDS$length,by = list(CDS$gene),sum),"CDS")
TE_combine = cbind(aggregate(TE$length,by=list(TE$gene),sum),"TE")
names(UTR_combine) <- names(exon_combine) <- names(CDS_combine) <- names(TE_combine)  <- c("gene","length","element")

#Calculate the intron length
intron_tem = merge(gene,exon_combine,by="gene")
intron_combine = data.frame(gene=intron_tem$gene,
                            length=intron_tem$length.x - intron_tem$length.y,
                            element="intron")

element_length_exon = merge(B73_core_epiallele,exon_combine,all.x = T)
element_length_intron = merge(B73_core_epiallele,intron_combine,all.x = T)
element_length_CDS = merge(B73_core_epiallele,CDS_combine,all.x = T)
element_length_UTR = merge(B73_core_epiallele,UTR_combine,all.x = T)
element_length_TE = merge(B73_core_epiallele,TE_combine,all.x = T)
element_length_TE$length[is.na(element_length_TE$length)] <- 0
element_length_UTR$length[is.na(element_length_UTR$length)] <-0
element_length_TE %>% group_by(epiallele) %>% summarise(mean(length))
element_length_UTR %>% group_by(epiallele) %>% summarise(sum(length>0))
element_length_combine = rbind(element_length_exon,element_length_intron,
                               element_length_CDS,element_length_UTR,
                               element_length_TE)
element_length_combine$length[is.na(element_length_combine$length)] = 0
element_length_combine$epiallele = factor(element_length_combine$epiallele,
                                          levels = c("UM","gbM","teM"))
element_length_combine$element = factor(element_length_combine$element,
                                          levels = c("UTR","TE","intron","CDS","exon"))
element_length_combine = element_length_combine[!is.na(element_length_combine$epiallele),]
element_length_combine = element_length_combine[!is.na(element_length_combine$element),]

#For TE frequency in intron
B73_TE = read.table("/Users/x/Desktop/Data/TE/B73_TE_intersect.bed")[,c(5,9)]
table(B73_TE$V9)
B73_TE = B73_TE[B73_TE$V9!="LTR",]
B73_TE = B73_TE[!duplicated(B73_TE),]
names(B73_TE) = c("gene","TE")

B73_TE_core = merge(B73_core_epiallele,B73_TE,by="gene")
df_B73_TE_core = table(B73_TE_core$epiallele,B73_TE_core$TE)/c(7541,710,12666)
df_B73_TE_core = data.frame(count=c(df_B73_TE_core),
                            epiallele=rep(row.names(df_B73_TE_core),9),
                            superfamily = rep(colnames(df_B73_TE_core),each=3)  )
df_B73_TE_core$epiallele = factor(df_B73_TE_core$epiallele,
                                  levels = c("UM","gbM","teM"))
df_B73_TE_core$superfamily = factor(df_B73_TE_core$superfamily,
                                    levels = c("CACTA","Harbinger","hAT","helitron","Mariner",
                                               "Copia","Gypsy","LINE","Mutator"))
#For UTR for diffferent epiallele
element_length_UTR$length[is.na(element_length_UTR$length)] = 0
UTR_count = table(element_length_UTR$length>0,element_length_UTR$epiallele)
df_UTR = data.frame(proportion = c(11940/(11940+726),7512/(7512+29),291/(419+291)),
                    methylation = c("UM","gbM","teM"))
df_UTR$methylation = factor(df_UTR$methylation,levels = c("UM","gbM","teM"))


#
CDS_TE = read.table("/Users/x/Desktop/Data/annotation/Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS_TE_intersect.bed")[,c(4,10)]
names(CDS_TE)<- c("gene","length")
CDS_TE_merge = aggregate(CDS_TE$length,by=list(CDS_TE$gene),sum)
names(CDS_TE_merge)<- c("gene","length")
df_CDS_TE_merge = merge(B73_core_epiallele,CDS_TE_merge,by="gene",all.x = T)
df_CDS_TE_merge$length[is.na(df_CDS_TE_merge$length)] = 0
freq = table(df_CDS_TE_merge$length>100,df_CDS_TE_merge$epiallele)
df_CDS_TE = data.frame(epiallele =c("gbM","teM","UM"),
                       proportion = c(396/sum(7145,396),215/sum(215,495),
                                      426/sum(426,12240)))
df_CDS_TE$epiallele = factor(df_CDS_TE$epiallele,levels = c("UM","gbM","teM"))
F2D=ggplot(df_CDS_TE,aes(x=epiallele,y=proportion)) + theme_bw() +
  geom_bar(stat = "identity") 

Hufford_data = read.table("/Users/x/Desktop/Data/expression/B73.tpm.txt",
                          header = T)
Hufford_data$status <- "tissue"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])>=1)==10] <- "expressed"
Hufford_data$status[rowSums(as.matrix(Hufford_data[,5:14])<1)==10] <- "silenced"
table(Hufford_data$status)
df_Hufford_data_core  = merge(B73_core_epiallele,Hufford_data,by="gene",all.x = T)
table(df_Hufford_data_core$epiallele,df_Hufford_data_core$status)
df_2E = data.frame(epiallele=rep(c("gbM","teM","UM"),3),
                   status = rep(c("constitutive","silent","tissue-specific"),each=3), 
                   count = c(table(df_Hufford_data_core$epiallele,df_Hufford_data_core$status)))

df_2E$epiallele = factor(df_2E$epiallele,levels = c("UM","gbM","teM"))
df_2E$status = factor(df_2E$status,levels = c("constitutive","silent","tissue-specific"))



df_TPM_core = data.frame(tissue=rep(names(df_Hufford_data_core)[6:15],each=dim(df_Hufford_data_core)[1]),
                         TPM = unlist(c(df_Hufford_data_core[,6:15])),
                         epiallele = rep(df_Hufford_data_core$epiallele,10))
df_TPM_core$epiallele = factor(df_TPM_core$epiallele,levels = c("UM","gbM","teM"))
df_TPM_core$tissue = factor(df_TPM_core$tissue,levels = names(df_Hufford_data_core)[c(11,12,13,9,8,10,15,14,6,7)])

#Comparison the average lengths of the UM/gbM/teM gene sets
element <- c("UM","gbM","teM")
context <- c("UTR","CDS","exon","intron","TE")
compare <- combn(element,2)
p_value <- c()
gr_context <- c()
gr_element <- c()
for (i in 1:3) {
  for (j in 1:5){
   p1 = wilcox.test(element_length$length[element_length$methylation==compare[1,i] & element_length$element== context[j]],
                element_length$length[element_length$methylation==compare[2,i] & element_length$element== context[j]])$p.value
   p_value <- c(p1,p_value)
   gr_context <- c(gr_context,paste0(compare[1,i],compare[2,i]))
   gr_element <- c(gr_element,context[j])
  }
}
cbind(p_value,gr_context,gr_element)



###########other length
UTR_length$length[is.na(UTR_length$length)] = 0
exon_length$length[is.na(exon_length$length)] = 0
intron_length$length[is.na(intron_length$length)] = 0
CDS_length$length[is.na(CDS_length$length)] = 0
TE_length$length[is.na(TE_length$length)] = 0

drop_na(UTR_length)%>% group_by(methylation)  %>% summarise(mean(length))
drop_na(exon_length) %>% group_by(methylation)  %>% summarise(mean(length))
drop_na(TE_length) %>% group_by(methylation)  %>% summarise(mean(length))
drop_na(CDS_length) %>% group_by(methylation)  %>% summarise(mean(length))
drop_na(intron_length) %>% group_by(methylation)  %>% summarise(mean(length))

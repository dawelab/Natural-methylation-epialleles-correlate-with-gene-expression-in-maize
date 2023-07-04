setwd("/Users/x/Desktop/Data/methylation/cgchgmtr")
#Read in the files containing epiallele information
#Read in the combine data frame 
path <- "/Users/x/Desktop/Data/combine/"
files <- list.files(path)
for (i in 1:26) {
  assign(
    gsub(".csv","",files[i]),
    read.csv(paste0(path,files[i]))
  )
}
NAM <- gsub(".all.csv","",files)

for (i in 1:26) {
  assign(
    paste0(NAM[i],".epi"),
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40) %>% 
      select(c("gene","class","epiallele")) %>% unique()
  )
}
#Loop the epiallele count into tidy data frame
epiallele_count <- data.frame(ambiguous = numeric(),
                              gbM = numeric(),   
                              teM = numeric(),
                              UM = numeric())
for (i in 1:26) {
  epiallele_count <- rbind(epiallele_count,table(select(get(paste0(NAM[i],".epi")),"epiallele")),
                           make.row.names = FALSE)
}
names(epiallele_count) <- c("ambiguous","gbM","teM","UM")
cbind(epiallele_count,NAM)
cbind(epiallele_count/rowSums(epiallele_count),NAM)
df_epi_count <- data.frame(count= c(unlist(epiallele_count, use.names = F)) ,
                           NAM = rep(NAM,4),
                           epiallele =rep(names(epiallele_count),each = 26))
df_epi_count$epiallele = factor(df_epi_count$epiallele,
                                   level = c("UM","gbM","teM","ambiguous"))
blank_group1 <- data.frame(epiallele = rep(c("UM","gbM","teM","ambiguous"),each=2),
                           x = 0 ,
                           y = c(0,1.5*10^4,0,10000,0,4000,0,6000))
#plotting
df_epi_count %>% ggplot(aes(x=NAM,y=count, fill = epiallele)) +
  geom_bar(stat = "identity") +
  geom_blank(data=blank_group1, aes(x=x,y=y)) +
  facet_wrap(~epiallele, nrow = 2 , scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("all gene")

##Subsetting the core gene
core_set <- paste0(NAM,".coreepi")
for (i in 1:26) {
  assign(
    core_set[i],
    get(gsub(".csv","",files[i])) %>% 
      filter(cCG >=40 & cCHG >= 40 & class == "Core Gene") %>% 
      select(c("gene","class","epiallele")) %>% unique()
  )
}

core_epiallele_count <- data.frame(ambiguous = numeric(),
                              gbM = numeric(),   
                              teM = numeric(),
                              UM = numeric())
for (i in 1:26) {
  temp <- get(core_set[i])   %>% select("epiallele") %>% table()
  core_epiallele_count <- rbind(core_epiallele_count,temp,
                           make.row.names = FALSE)
}
names(core_epiallele_count) <- c("ambiguous","gbM","teM","UM")
df_epi_count_core <- data.frame(count= c(unlist(core_epiallele_count, use.names = F)) ,
                           NAM = rep(NAM,4),
                           epiallele =rep(names(core_epiallele_count),each = 26))
df_epi_count_core$epiallele = factor(df_epi_count$epiallele,
                                     level = c("UM","gbM","teM","ambiguous"))

blank_group2 <- data.frame(epiallele = rep(c("UM","gbM","teM","ambiguous"),each=2),
                           x = 0 ,
                           y = c(0,1.5*10^4,0,10000,0,4000,0,6000))

#plotting
df_epi_count_core %>% ggplot(aes(x=NAM,y=count, fill = epiallele)) +
  geom_bar(stat = "identity") +
  geom_blank(data=blank_group2,aes(x=x,y=y)) +
  facet_wrap(~epiallele, nrow = 2 , scales = "free_y") +
  theme_bw() + theme(axis.text.x  = element_text(angle = 90)) +
  ggtitle("core gene")

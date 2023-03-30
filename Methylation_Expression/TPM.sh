cd /scratch/yz77862/RNA/self/merge
ls *bam | sed 's/.bam//g' > list
ml R
list=/scratch/yz77862/RNA/self/counts/list
while read LIB_LINE;do
   GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f1)
   OUT=/scratch/yz77862/shell/tpm/${GENOME}_self.R
    echo '#!/usr/bin/env Rscript' >> ${OUT}  
    echo "setwd(\"/scratch/yz77862/RNA/self/tpm\")"  >> ${OUT} 
    echo "data=read.table(\"/scratch/yz77862/RNA/self/counts/${GENOME}\",header=T)"  >> ${OUT}  
    echo "#To calculate the tpm "  >> ${OUT} 
    echo " rpkm <- function(counts, lengths) {"  >> ${OUT} 
    echo " rate <- counts / lengths"  >> ${OUT} 
    echo "rate / sum(counts) * 1e9 "  >> ${OUT} 
    echo " }"  >> ${OUT} 
    echo " tpm <- function(counts, lengths) {"  >> ${OUT} 
    echo " rate <- counts / lengths"  >> ${OUT} 
    echo " rate / sum(rate) * 1e6"  >> ${OUT} 
    echo " }"  >> ${OUT} 
    echo "data\$tpm=tpm(data[,7],data[,6])"  >> ${OUT} 
    echo "data\$rpkm=rpkm(data[,7],data[,6])"  >> ${OUT} 
    echo "df_tpm=data.frame(data[,c(1,8)])"  >> ${OUT} 
    echo "df_rpkm=data.frame(data[,c(1,9)])"  >> ${OUT} 
    echo "write.table(df_tpm, file=\"${GENOME}_tpm.txt\", sep=\"\\t\", row.names=FALSE, quote=FALSE)"  >> ${OUT} 
    echo "write.table(df_rpkm, file=\"${GENOME}_rpkm.txt\", sep=\"\\t\", row.names=FALSE, quote=FALSE)"  >> ${OUT} 
    chmod +x ${OUT}
    
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)

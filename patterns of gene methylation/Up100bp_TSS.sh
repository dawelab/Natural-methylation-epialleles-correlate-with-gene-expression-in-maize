#To make the bed files for cgmaptools mtr command
cd /scratch/yz77862/gene_methylation/mtr/gene
for methy in UM gbM teM;do
    awk '$10=="'${methy}'"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt | awk '$4=="+"' | awk '{print $1,$2-101,$2-1}' OFS="\t" > B73.${methy}_up.bed
    awk '$10=="'${methy}'"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt | awk '$4=="1"' | awk '{print $1,$3+1,$3+101}' OFS="\t" > B73.${methy}_down.bed
    cat B73.${methy}_up.bed B73.${methy}_down.bed | sort -b -k1,1 -k2,2n -k3,3n > B73.${methy}_TSS.bed
    rm B73.${methy}_up.bed B73.${methy}_down.bed 
done 
    
          
#Use the cgmaptools mtr command to calculate methylation levels
ml CGmapTools/0.1.2-foss-2019b

input=/scratch/yz77862/gene_methylation/CGmap/merge/B73_B73v5_merge.CGmap
region=/scratch/yz77862/gene_methylation/mtr/gene
out=/scratch/yz77862/gene_methylation/mtr/upTSS

for methy in UM gbM teM;do
    for context in CG CHG CHH;do
        awk '$4=='"${context}"'' ${input}| cgmaptools mtr -r ${region}/B73.${methy}_TSS.bed >${out}/${GENOME}.${methy}.${context}.TSS.100.mtr
done
    done

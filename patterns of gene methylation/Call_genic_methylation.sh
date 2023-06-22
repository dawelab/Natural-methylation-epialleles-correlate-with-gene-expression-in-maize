#Calculate the gene methylation level on gene region

merged_list=/scratch/yz77862/gene_methylation/gene_annotation/NAM.list
mkdir -p /scratch/yz77862/gene_methylation/CGmap/gene

while read LIB_LINE; do
    Founder=$(echo ${LIB_LINE} | cut -d ' ' -f1)
    GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f3)
   
    CGmap=/scratch/yz77862/gene_methylation/CGmap/${Founder}
    CDS=/scratch/yz77862/gene_methylation/gene_annotation/${GENOME}.CDS
    outpath=/scratch/yz77862/gene_methylation/CGmap/gene/

    zcat  ${CGmap}| cgmaptools select region -r ${CDS}  > ${outpath}/${GENOME}.CDS.CGmap

done < <(cut -f1,2,3,4 ${merged_list} | grep -v 'skip' | sort -u)

cgmaptools mtr -i WG.CGmap.gz -r region.bed -o WG.mtr.gz

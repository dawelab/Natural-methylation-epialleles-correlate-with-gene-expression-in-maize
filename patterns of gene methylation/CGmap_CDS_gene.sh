#convert gff3 file to bed file
list=/scratch/yz77862/gene_methylation/gene_annotation/canonical/list
while read GENOME;do
gff=/scratch/yz77862/gene_methylation/gene_annotation/canonical/${GENOME}.gff3
sed 's/chr//g' ${gff}| awk '$3=="gene"' | awk '{print $1,$4-1,$5,$7,$9}' OFS="\t" | sed 's/;biotype=protein_coding;logic_name=cshl_gene//g' | sed 's/ID=//g' | sort -b -k1,1 -k2,2n -k3,3n > ${GENOME}.gene.bed
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)

#convert gft file to bed file with attached geneID
while read GENOME;do
awk '{print $1,$4,$5,$12}' OFS="\t" /scratch/yz77862/gene_methylation/gene_annotation/canonical/${GENOME}.gtf | sed 's/"//g' | sed 's/;//g' > ${GENOME}.gene.bed
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)


#Make the canonical CDS bed file
while read GENOME;do
sed 's/chr//g' ${GENOME}.gff3 | awk '$3=="CDS"'  | awk '{print $1,$4-1,$5}' OFS="\t" > ${GENOME}.CDS.bed
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)

#Use cgmaptools select region to select CDS region from methylomes
list=/scratch/yz77862/gene_methylation/gene_annotation/canonical/list3
while read LIB_LINE; do
   Founder=$(echo ${LIB_LINE} | cut -d ' ' -f1)
   GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f2)
   OUT=/scratch/yz77862/shell/${Founder}.sh
  echo '#!/bin/bash' >> ${OUT} 
  echo "#SBATCH --job-name=${Founder}" >> ${OUT} 
  echo "#SBATCH --partition=batch" >> ${OUT} 
  echo "#SBATCH --ntasks=1" >> ${OUT} 
  echo "#SBATCH --mem=180G" >> ${OUT} 
  echo "#SBATCH --time=8:00:00" >> ${OUT} 
  echo "#SBATCH --export=NONE" >> ${OUT} 
  echo "#SBATCH --output=${Founder}_Index_star.out" >> ${OUT} 
  echo "#SBATCH --error=${Founder}_Index_star.err" >> ${OUT}
  echo "#SBATCH --error=${Founder}_Index_star.err" >> ${OUT}
  echo "ml CGmapTools/0.1.2-foss-2019b" >> ${OUT}
  echo "cgmaptools select region -i /scratch/yz77862/gene_methylation/CGmap/merge/${GENOME} -r /scratch/yz77862/gene_methylation/gene_annotation/canonical/${Founder}.CDS.bed2 > /scratch/yz77862/gene_methylation/CGmap/CDS/${Founder}.CDS.CGmap" >> ${OUT}
  echo "cd /scratch/yz77862/gene_methylation/mtr/gene" >> ${OUT}
  echo "awk '\$4==\"CG\"' /scratch/yz77862/gene_methylation/CGmap/CDS/${Founder}.CDS.CGmap| cgmaptools mtr -r /scratch/yz77862/gene_methylation/gene_annotation/canonical/${Founder}.gene.bed -o ${Founder}.gene.CG.mtr" >> ${OUT}
  echo "awk '\$4==\"CHG\"' /scratch/yz77862/gene_methylation/CGmap/CDS/${Founder}.CDS.CGmap| cgmaptools mtr -r /scratch/yz77862/gene_methylation/gene_annotation/canonical/${Founder}.gene.bed -o ${Founder}.gene.CHG.mtr" >> ${OUT}

done < <(cut -f1,2 ${list} | grep -v 'skip' | sort -u)

list=/scratch/yz77862/gene_methylation/gene_annotation/canonical/list3
while read LIB_LINE; do
   Founder=$(echo ${LIB_LINE} | cut -d ' ' -f1)
   GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f2)
   paste ${Founder}.gene.CG.mtr ${Founder}.gene.CHG.mtr ${Founder}.gene.bed| awk '{print $1,$2,$3,$18,$6,$5,$13,$12,$19}' OFS="\t" > ${Founder}.gene.cgchg.mtr
   done < <(cut -f1,2 ${list} | grep -v 'skip' | sort -u)

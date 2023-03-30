#!/bin/bash
#To create the index file for NAM founder
#Authors: Yibing Zeng
#Date: 02/2/2022

merged_list=/scratch/yz77862/genome/list

while read LIB_LINE; do
   Founder=$(echo ${LIB_LINE} | cut -d ' ' -f1)
   GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f2)

  OUT=/scratch/yz77862/shell/star_index/${Founder}.sh
  echo '#!/bin/bash' >> ${OUT} 
  echo "#SBATCH --job-name=${Founder}" >> ${OUT} 
  echo "#SBATCH --partition=highmem_p" >> ${OUT} 
  echo "#SBATCH --ntasks=1" >> ${OUT} 
  echo "#SBATCH --mem=180G" >> ${OUT} 
  echo "#SBATCH --time=8:00:00" >> ${OUT} 
  echo "#SBATCH --export=NONE" >> ${OUT} 
  echo "#SBATCH --output=${Founder}_Index_star.out" >> ${OUT} 
  echo "#SBATCH --error=${Founder}_Index_star.err" >> ${OUT}
  echo " " >> ${OUT} 
  echo "ml STAR/2.7.2b-GCC-8.3.0" >> ${OUT} 
  echo "ml Cufflinks/2.2.1-foss-2019b" >> ${OUT} 
  echo " " >> ${OUT} 
  echo "#To make the variable for the gtf files" >> ${OUT} 
  echo "gff=/scratch/yz77862/gene_methylation/gene_annotation/${GENOME}.gff3" >> ${OUT}  
  echo "sed 's/chr//g' \${gff} | gffread -T -o /scratch/yz77862/gene_methylation/gene_annotation/${GENOME}.gtf" >> ${OUT}  
  echo "gtf=/scratch/yz77862/gene_methylation/gene_annotation/${GENOME}.gtf" >> ${OUT} 
  echo ""  >> ${OUT} 
  echo "#The pathway to fasta file"  >> ${OUT} 
  echo "fasta=/scratch/yz77862/genome/${Founder}.fasta" >> ${OUT} 
  echo " "  >> ${OUT} 
  echo "#The pathway to store the index"  >> ${OUT} 
  echo "mkdir -p /scratch/yz77862/index/${Founder}/STAR"  >> ${OUT} 
  echo "dir=/scratch/yz77862/index/${Founder}/STAR" >> ${OUT} 
  echo " "  >> ${OUT} 
  echo "STAR \\"  >> ${OUT} 
  echo "     --runThreadN 8 \\"  >> ${OUT} 
  echo "     --runMode genomeGenerate \\"  >> ${OUT} 
  echo "     --genomeDir \${dir} \\"  >> ${OUT} 
  echo "     --genomeFastaFiles \${fasta} \\"  >> ${OUT} 
  echo "     --sjdbGTFfile \${gtf}"  >> ${OUT} 
sbatch ${OUT}
done < <(cut -f1,2 ${merged_list} | grep -v 'skip' | sort -u)

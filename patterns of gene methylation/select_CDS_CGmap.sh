#To select CDS region

cd /scratch/yz77862/gene_methylation/CGmap/CDS
genome_list=/scratch/yz77862/gene_methylation/gene_annotation/NAM.list
while read LIB_LINE; do
   Founder=$(echo ${LIB_LINE} | cut -d ' ' -f1)
   genome=$(echo ${LIB_LINE} | cut -d ' ' -f3)

  OUT=/scratch/yz77862/gene_methylation/CGmap/CDS/${genome}.extract_CDS.sh

  echo '#!/bin/bash' >> ${OUT}
  echo "#SBATCH --job-name=TE2gff" >> ${OUT}
  echo "#SBATCH --partition=batch" >> ${OUT}
  echo "#SBATCH --ntasks=1" >> ${OUT}
  echo "#SBATCH --mem=100G" >> ${OUT}
  echo "#SBATCH --time=00:03:00" >> ${OUT}
  echo "#SBATCH --export=NONE" >> ${OUT}
  echo "#SBATCH --output=${genome}.out" >> ${OUT}
  echo "#SBATCH --error=${genome}.err" >> ${OUT}
  echo "" >> ${OUT}
  echo "" >> ${OUT}
  echo "ml CGmapTools/0.1.2-foss-2019b" >> ${OUT}
  echo "cd /scratch/yz77862/gene_methylation/CGmap/CDS" >> ${OUT}
  echo "CGmap=/scratch/yz77862/gene_methylation/CGmap/${Founder}" >> ${OUT}
  echo "CDS=/scratch/yz77862/gene_methylation/gene_annotation/${genome}.CDS" >> ${OUT}
  echo "outpath=/scratch/yz77862/gene_methylation/CGmap/CDS/" >> ${OUT}
  echo "" >> ${OUT}
  echo "zcat  \${CGmap}| cgmaptools select region -r \${CDS}  > \${outpath}/${genome}.CDS.CGmap" >> ${OUT}

done < <(cut -f1,2,3,4 ${genome_list} | grep -v 'skip' | sort -u)

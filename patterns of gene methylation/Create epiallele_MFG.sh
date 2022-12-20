list=/scratch/yz77862/gene_methylation/mtr/gene/MFG.list #The list of region file from preMFG.sh 
while read Region; do
    OUT=/scratch/yz77862/gene_methylation/mfg/shell/${Region}.sh
    echo '#!/bin/bash' >> ${OUT}
    echo "#SBATCH --job-name=${Region}.MFG" >> ${OUT}
    echo "#SBATCH --partition=batch" >> ${OUT}
    echo "#SBATCH --ntasks=1" >> ${OUT}
    echo "#SBATCH --mem=100G" >> ${OUT}
    echo "#SBATCH --time=030:00:00" >> ${OUT}
    echo "#SBATCH --export=NONE" >> ${OUT}
    echo "#SBATCH --output=${Region}.out" >> ${OUT}
    echo "#SBATCH --error=${Region}.err" >> ${OUT}
    echo "" >> ${OUT}
    echo "ml CGmapTools/0.1.2-foss-2019b" >> ${OUT}
    echo "cd /scratch/yz77862/gene_methylation/mtr/gene" >> ${OUT}
    echo "I=/scratch/yz77862/gene_methylation/CGmap/merge/B73_B73v5_merge.CGmap" >> ${OUT}
    echo "R=${Region}_sort.bed" >> ${OUT}
    echo "out=/scratch/yz77862/gene_methylation/mfg/context" >> ${OUT}
    echo "    for context in CG CHG CHH;do" >> ${OUT}
    echo "" >> ${OUT}
    echo "        cgmaptools mfg -i \${I} -r \${R} -c 1 -x \${context}> \${out}/${Region}.\${context}.mfg" >> ${OUT}
    echo "" >> ${OUT}    
    echo "done" >> ${OUT}
        
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)

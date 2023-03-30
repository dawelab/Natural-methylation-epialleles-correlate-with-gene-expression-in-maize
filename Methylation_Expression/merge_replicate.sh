list=/project/krdlab/Yibing/RNA_map/fastq/bam_merge.list
while read INPUT; do
OUT=/scratch/yz77862/shell/merge_bam/${INPUT}_merge_bam_self.sh
    echo '#!/bin/bash' >> ${OUT}
    echo "#SBATCH --job-name=${INPUT}_merge_bam"  >> ${OUT}
    echo "#SBATCH --partition=batch"  >> ${OUT}
    echo "#SBATCH --ntasks=1"  >> ${OUT}
    echo "#SBATCH --mem=100G"  >> ${OUT}
    echo "#SBATCH --time=5:00:00"  >> ${OUT}
    echo "#SBATCH --export=NONE"  >> ${OUT}
    echo "#SBATCH --output=${INPUT}_merge_bam.out"  >> ${OUT}
    echo "#SBATCH --error=${INPUT}_merge_bam.err"  >> ${OUT}
    echo " "  >> ${OUT}
    echo "ml SAMtools/0.1.19-foss-2019b"  >> ${OUT}
    echo "cd /scratch/yz77862/RNA/self/merge"  >> ${OUT}

    echo "samtools merge ${INPUT}.merged.bam /scratch/yz77862/RNA/self/round2/${INPUT}_*_round-2Aligned.sortedByCoord.out.bam"  >> ${OUT}
done < <(cut -f1 ${list} | grep -v 'skip' | sort -u)

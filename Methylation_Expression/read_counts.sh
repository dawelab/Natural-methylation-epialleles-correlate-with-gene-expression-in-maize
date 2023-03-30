list=/scratch/yz77862/shell/tpm/list6
while read LIB_LINE;do
   GTF=$(echo ${LIB_LINE} | cut -d ' ' -f2)
   GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f1)

   OUT=/scratch/yz77862/shell/tpm/${GENOME}_count_self.sh
    echo '#!/bin/bash'   >> ${OUT}  
    echo "#SBATCH --job-name=${GENOME}_self.count"   >> ${OUT} 
    echo "#SBATCH --partition=batch"  >> ${OUT}  
    echo "#SBATCH --nodes=1"  >> ${OUT}  
    echo "#SBATCH --ntasks=1"  >> ${OUT}  
    echo "#SBATCH --cpus-per-task=1"  >> ${OUT}  
    echo "#SBATCH --mem=45G"  >> ${OUT}  
    echo "#SBATCH --time=00:30:00"  >> ${OUT}  
    echo "#SBATCH --output=${GENOME}_self.out"  >> ${OUT}  
    echo "#SBATCH --error=${GENOME}_self.err"  >> ${OUT}  
    echo "cd /scratch/yz77862/RNA/self/merge/"  >> ${OUT} 
    echo " "  >> ${OUT} 
    echo "Bam=${GENOME}"  >> ${OUT}
    echo "gtf=/scratch/yz77862/gene_methylation/gene_annotation/${GTF}"  >> ${OUT} 
    echo " "  >> ${OUT}  
    echo "/home/yz77862/apps/subread-1.6.0-Linux-x86_64/bin/featureCounts  -a \${gtf} -o /scratch/yz77862/RNA/self/counts/${GENOME}.counts \${Bam}" >> ${OUT} 
    sbatch ${OUT}
done < <(cut -f1,2,3 ${list} | grep -v 'skip' | sort -u)

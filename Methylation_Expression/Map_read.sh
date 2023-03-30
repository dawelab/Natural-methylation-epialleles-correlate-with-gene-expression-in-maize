fq_list=/project/krdlab/Yibing/RNA_map/fastq/self_map.list
while read INPUT; do
   Founder=$(echo ${INPUT} | cut -d ' ' -f2)
   GENOME=$(echo ${INPUT} | cut -d ' ' -f1)

OUT=/scratch/yz77862/shell/fastq_bam/${GENOME}_self.sh
    echo '#!/bin/bash'  >> ${OUT} 
    echo "#SBATCH --job-name=${GENOME}_self_mapping"   >> ${OUT}            
    echo "#SBATCH --partition=batch"   >> ${OUT} 
    echo "#SBATCH --nodes=1"   >> ${OUT}                  
    echo "#SBATCH --ntasks=1"   >> ${OUT}              
    echo "#SBATCH --cpus-per-task=18"   >> ${OUT}          
    echo "#SBATCH --mem=80G"   >> ${OUT}                  
    echo "#SBATCH --time=010:00:00"   >> ${OUT}             
    echo "#SBATCH --output=${GENOME}_fq_bam.out"   >> ${OUT}         
    echo "#SBATCH --error=${GENOME}_fq_bam.err"   >> ${OUT}         
    echo "#SBATCH --mail-user=yz77862@uga.edu"   >> ${OUT}  
    echo " "  >> ${OUT}  
    echo "ml STAR/2.7.2b-GCC-8.3.0" >> ${OUT}  
    echo " "  >> ${OUT}
    echo "mkdir -p  /scratch/yz77862/RNA/self/round1"  >> ${OUT}
    echo "cd /scratch/yz77862/RNA/self/round1" >> ${OUT}
    echo " "  >> ${OUT}
    echo "thread=18"  >> ${OUT}  
    echo "index=/scratch/yz77862/index/${Founder}/STAR/"  >> ${OUT}  
    echo "read1=/scratch/yz77862/RNA/fastq/${GENOME}_R1.fq.gz"  >> ${OUT}  
    echo "read2=/scratch/yz77862/RNA/fastq/${GENOME}_R2.fq.gz"  >> ${OUT}  
    echo " "  >> ${OUT}  
    echo "STAR \\"  >> ${OUT}    
    echo "--runMode alignReads \\"  >> ${OUT}  
    echo "--twopassMode Basic  \\"  >> ${OUT}  
    echo "​--runThreadN \$thread \\"  >> ${OUT}  
    echo "--genomeDir \${index}  \\"  >> ${OUT}  
    echo " --readFilesIn \${read1} \${read2} \\"  >> ${OUT}  
    echo "--readFilesCommand zcat \\"  >> ${OUT}  
    echo "--outSAMtype None \\"  >> ${OUT}  
    echo "--outFileNamePrefix ${GENOME}"  >> ${OUT}  
    echo " "  >> ${OUT}
    echo "mkdir -p  /scratch/yz77862/RNA/self/round2"  >> ${OUT}
    echo "cd /scratch/yz77862/RNA/self/round2"  >> ${OUT}
    echo " " >> ${OUT}
    echo "SJ=/scratch/yz77862/RNA/self/round1/${GENOME}_STARpass1/SJ.out.tab"  >> ${OUT}
    echo "STAR \\"  >> ${OUT}
    echo "--genomeDir \${index} \\"  >> ${OUT}
    echo "--runThreadN \${thread} \\"  >> ${OUT}
    echo "--sjdbFileChrStartEnd \${SJ} \\"  >> ${OUT}
    echo "--runMode alignReads \\"  >> ${OUT}
    echo "--readFilesIn \${read1} \${read2} \\"  >> ${OUT}
    echo "--readFilesCommand zcat \\"  >> ${OUT}
    echo "--outSAMattributes All \\"  >> ${OUT}
    echo "--outSAMmapqUnique 10 \\"  >> ${OUT}
    echo "--outFilterMismatchNmax 0 \\"  >> ${OUT}
    echo "--outFileNamePrefix ${GENOME}_round-2 \\"  >> ${OUT}
    echo "--outBAMsortingThreadN 4 \\"  >> ${OUT}
    echo "--limitBAMsortRAM 5594394835 \\"  >> ${OUT}
    echo "--outSAMtype BAM SortedByCoordinate \\"  >> ${OUT}
    echo "--outWigType bedGraph read1_5p"  >> ${OUT}
#sbatch ${OUT}
done < <(cut -f1,2 ${fq_list} | grep -v 'skip' | sort -u)

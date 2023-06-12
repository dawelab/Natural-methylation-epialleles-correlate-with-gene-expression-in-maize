while read GENOME;do
    echo '#!/bin/bash' >> ${OUT}
    echo "#SBATCH --job-name=extract_CDS" >> ${OUT}
    echo "#SBATCH --partition=batch" >> ${OUT}
    echo "#SBATCH --ntasks=1" >> ${OUT}
    echo "#SBATCH --mem=100G" >> ${OUT}
    echo "#SBATCH --time=030:00:00" >> ${OUT}
    echo "#SBATCH --export=NONE" >> ${OUT}
    echo "#SBATCH --output=${GENOME}.out" >> ${OUT}
    echo "#SBATCH --error=${GENOME}.err" >> ${OUT}
    echo "" >> ${OUT}
    echo "" >> ${OUT}
    echo "ml CGmapTools/0.1.2-foss-2019b" >> ${OUT}
    echo "cd " >> ${OUT}
    echo "input=" >> ${OUT}
    echo "region=" >> ${OUT}
    echo "out=" >> ${OUT}
    echo "" >> ${OUT}
    echo "cgmaptools select region -i \${input} -r \${region} > \${out}/${GENOME}.CDS.CGmap" >> ${OUT}
        
done < <(cut -f1 ${list})   

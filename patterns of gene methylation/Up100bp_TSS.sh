   #!/bin/bash
#SBATCH --job-name=B73up100TSS
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --export=NONE 
#SBATCH --output=B73up100TSS.out
#SBATCH --error=B73up100TSS.err
    
ml CGmapTools/0.1.2-foss-2019b
cd /scratch/yz77862/gene_methylation
#B73gff3=Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3
#make the bed file for 100bp upstream TSS
#cut -f1  ${B73gff3} -d ';'| awk '$3=="gene"' | sed 's/ID=//g' | sed 's/chr//g' | awk '{if ($7=="+") {print $1,$4-101,$4-1,$7,$9} else if ($7=="-") {print $1,$5+1,$5+100,$7,$9}}' OFS="\t" > Zm-B73-REFERENCE-NAM-5.0.1.canon.100TSS.bed
B73up100bed=Zm-B73-REFERENCE-NAM-5.0.1.canon.100TSS.bed
CGmap=B73_B73v5_merge.CGmap
awk '$4=="CG"' ${CGmap} | cgmaptools mtr -r ${B73up100bed} > B73.canno.100upTSS.CG.mtr
awk '$4=="CHG"' ${CGmap} | cgmaptools mtr -r ${B73up100bed} > B73.canno.100upTSS.CHG.mtr
awk '$4=="CHH"' ${CGmap} | cgmaptools mtr -r ${B73up100bed} > B73.canno.100upTSS.CHH.mtr

paste Zm-B73-REFERENCE-NAM-5.0.1.canon.100TSS.bed B73.canno.100upTSS.CG.mtr | awk '{print $5,$11,"CG"}' > B73.canno.100upTSS.mCG
paste Zm-B73-REFERENCE-NAM-5.0.1.canon.100TSS.bed B73.canno.100upTSS.CHG.mtr | awk '{print $5,$11,"CHG"}' > B73.canno.100upTSS.mCHG
paste Zm-B73-REFERENCE-NAM-5.0.1.canon.100TSS.bed B73.canno.100upTSS.CHH.mtr | awk '{print $5,$11,"CHH"}' > B73.canno.100upTSS.mCHH

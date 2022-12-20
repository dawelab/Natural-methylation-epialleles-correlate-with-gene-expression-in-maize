#Prepare files for CGmaptools MFG analysis
#!/bin/bash

ml  CGmapTools/0.1.2-foss-2019b
I=/scratch/yz77862/gene_methylation/mtr/gene/  #Put in the the bed (anotation) files
#To generate the UM bed file, I used the one that generated in epiallle type:
#awk '$10=="UM"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.gene.mtr.ID.type.txt | awk '{print $1,$2,$3,$4}' OFS="\t" > B73.UM.bed
prefix= #Put in the name of you want for the output. Name the file!
let withingene=2000  #Put in the number about how deep you want to go into the gene
let bin=100    #Put in the number of bin size
let flank=3000  #Put in the number of the flanking region
let n=(${withingene}+${flank})/$bin
awk '$4=="+"' ${I}| awk '{print $1,$2-'$flank'-1,$2+'$withingene'-1,"+"}' OFS="\t" > ${prefix}_up1.bed
awk '$4=="-"' ${I} |  awk '{print $1,$3-'$withingene',$3+'$flank',"-"}' OFS="\t" > ${prefix}_down1.bed
cat  ${prefix}_up1.bed ${prefix}_down1.bed | sort -b -k1,1 -k2,2n -k3,3n> ${prefix}.UP.bed
cgmaptools bed2fragreg -i ${prefix}.UP.bed -n $n -o ${prefix}_UP_final.bed
sort -b -k1,1 -k2,2n -k3,3n ${prefix}_UP_final.bed > ${prefix}_UP_sort.bed
rm ${prefix}.UP.bed ${prefix}_UP_final.bed ${prefix}_up1.bed ${prefix}_down1.bed
#Calculating the downstream
awk '$4=="+"' ${I} | awk '{print $1,$3-'$withingene',$3+'$flank',"+"}' OFS="\t" > ${prefix}_up2.bed
awk '$4=="-"' ${I} |  awk '{print $1,$2-'$flank'-1,$2+'$withingene'-1,"-"}' OFS="\t" > ${prefix}_down2.bed
cat ${prefix}_up2.bed ${prefix}_down2.bed | sort -b -k1,1 -k2,2n -k3,3n > ${prefix}_DOWN.bed
cgmaptools bed2fragreg -i ${prefix}_DOWN.bed -n $n -o ${prefix}_DOWN_final.bed
sort -b -k1,1 -k2,2n -k3,3n ${prefix}_DOWN_final.bed > ${prefix}_DOWN_sort.bed
rm ${prefix}_up2.bed  ${prefix}_down2.bed ${prefix}_DOWN_final.bed ${prefix}_DOWN.bed ${prefix}.bed

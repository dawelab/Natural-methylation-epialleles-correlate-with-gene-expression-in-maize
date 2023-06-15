#!/usr/bin/env bash

#This is how to define the methylation type for everygenes
#It requires the file with CG/CHG methylation levels and their #C

#The list of NAM founder used in this analysis
genome_list=/Users/x/Desktop/Data/methylation/cgchgmtr/list

#The organized input files as columned as chromosome number, start, end, strand, coverage for CG, methylation level for CG, coverage for CHG, methylation level for CHG:
#1       41823   47427   +       0.2950  172     0.0008  166     Zm00036ab000010
#1       48390   53991   -       0.3358  104     0.0023  166     Zm00036ab000020
#1       118413  118791  +       0.4074  8       0.0000  35      Zm00036ab000050
#1       197942  199261  -       0.7951  65      0.6299  74      Zm00036ab000070
#1       199291  206669  -       0.8120  78      0.6098  75      Zm00036ab000080
#1       208099  212734  -       0.6560  30      0.0000  48      Zm00036ab000090
#1       212904  213736  -       0.0000  3       NA      0       Zm00036ab000100


while read LIB_LINE; do
   genome=$(echo ${LIB_LINE} | cut -d ' ' -f3)
   #Defining all the parameters in determining the gene methylation type
   
   cCG=40         # The thredholds of the number for the cytosines in CDS regions at CG context
   cCHG=40        # The thredholds of the number for the cytosines in CDS regions at CHG context
   mCG1=0.05       #The thredholds of the unmethylated region at CG context
   mCHG1=0.05      #The thredholds of the unmethylated region at CHG context
   mCG2=0.2      #The thredholds of the gene body methylated region at CG context
   mCHG2=0.05     #The thredholds of the gene body methylated region at CHG context
   mCG3=0.4      #The CG methylation level thredholds of teM
   mCHG3=0.4     #The CHG methylation level thredholds of teM
   
   #In defining the unmethylated genes, it requires the high coverages for the genes and low methylation levels at CG and CHG context (<=0.2)  
   #In defining the gene body methylated genes, it requires the high coverages for the genes and low methylation levels CHG context (<=0.2) but high methylation level atat CG context (>=0.4) 
   #In defining the transposable like genes, it requires the high coverages for the genes and high methylation levels at both CG and CHG context (>=0.4)  

   awk '{if (($6>='"${cCG}"') && ($8>='"${cCHG}"') && ($5<='"${mCG1}"') && ($7<='"${mCHG1}"')) {$10=""; print $0,"UM"} \
   else if (($6>='"${cCG}"') && ($8>='"${cCHG}"') && ($5>='"${mCG2}"') && ($7<='"${mCHG1}"')) {$10=""; print $0,"gbM"} \
   else if (($6>='"${cCG}"') && ($8>='"${cCHG}"') && ($5>='"${mCG3}"') && ($7>='"${mCHG3}"')) {$10=""; print $0,"teM"} \
   else if (($6>='"${cCG}"') && ($8>='"${cCHG}"') && ($5>'"${mCG1}"') && ($5<'"${mCG2}"') && ($7<'"${mCHG1}"')) {$10=""; print $0,"interCG"} \
   else if (($6>='"${cCG}"') && ($8>='"${cCHG}"') && ($5>'"${mCG1}"') && ($5<'"${mCG3}"') && ($7>'"${mCHG1}"') && ($7<'"${mCHG3}"')) {$10=""; print $0,"interCHG"} \
   else {$10=""; print $0,"ambiguous"}}' OFS="\t" ${genome}.gene.gene.mtr.ID.type.txt > ${genome}.gene.mtr.ID.type
   
 
done < <(cut -f1,2,3,4 ${genome_list} | grep -v 'skip' | sort -u)

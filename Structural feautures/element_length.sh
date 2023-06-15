#All of the output data will be stored as tidy dataframe as: element_type length geneID
##Generate a gff3 file with intron information using genome tools
ml  GenomeTools/1.6.1-GCC-10.2.0
gt gff3 -addintrons -retainids -sortlines Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 > Zm-B73-REFERENCE-NAM-5.0.1.canon.addintron.gff3

#use the canonical B73 annotation file as the input to calculate the element length of CDS, exon and UTR
awk '$3=="exon" || $3=="intron" ||$3 =="CDS" || $3=="five_prime_UTR" || $3=="three_prime_UTR"' Zm-B73-REFERENCE-NAM-5.0.1.canon.addintron.gff3 |\
cut -f1 -d ';' | sed 's/ID=//g' | sed 's/Parent=//g' | sed 's/_.00.//g' | \
awk '{print $3,$5-$4+1,$9}' OFS="\t" > B73_CDS_exon_intron_UTR.txt


#generate bed file of intron feature from gff3 file
awk '$3=="intron"' Zm-B73-REFERENCE-NAM-5.0.1.canon.addintron.gff3 | cut -f1 -d ';' | \
sed 's/chr//g' | sed 's/Parent=//g' | sed 's/_.00.//g' |\
awk '{print $1,$4-1,$5,".",".",$7,$9}' OFS="\t" | awk '$3!="."' | \
sort -b -k1,1 -k2,2n -k3,3n  > Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.bed

ml BEDTools
#use the TE annotation file as the input to generate the TE bed file
#Merge the TE feautures that ovelapped with each other
awk '{print $1,$4-1,$5,$3}' OFS="\t" Zm-B73-REFERENCE-NAM-5.0.TE.gff3 | sed 's/chr//g' | sort -b -k1,1 -k2,2n -k3,3n > Zm-B73-REFERENCE-NAM-5.0.TE.bed
bedtools merge -i Zm-B73-REFERENCE-NAM-5.0.TE.bed > Zm-B73-REFERENCE-NAM-5.0.TE.merged.bed

#Intersect the intron bed file and TE merged bed file
bedtools intersect -wo -a Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.bed -b Zm-B73-REFERENCE-NAM-5.0.TE.merged.bed > Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.bed

#Calculate the TE length intersect with intron
awk '{print "TE",$10-$9,$7}' OFS="\t" Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.bed > Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.txt

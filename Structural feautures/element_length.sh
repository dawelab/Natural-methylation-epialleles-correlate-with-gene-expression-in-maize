#All of the output data will be stored as tidy dataframe as: element_type length geneID
#use the canonical B73 annotation file as the input to calculate the element length of CDS, exon and UTR
awk '$3=="exon" || $3 =="CDS" || $3=="five_prime_UTR" || $3=="three_prime_UTR"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 |\
cut -f1 -d ';' | sed 's/ID=//g' | sed 's/Parent=//g' | sed 's/_.00.//g' | \
awk '{print $3,$5-$4+1,$9}' OFS="\t" > B73_CDS_exon_UTR.txt

#Generate the intron file
ml BEDTools/2.30.0-GCC-10.2.0
#generate bed file of gene feature from gff3 file
awk '$3=="gene"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 | cut -f1 -d ';' | \
sed 's/ID=//g' | sed 's/chr//g' | \
awk '{print $1,$4-1,$5,".",".",$7,$9}' OFS="\t" | \
sort -b -k1,1 -k2,2n -k3,3n  > Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.bed

#generate bed file of exon feature from gff3 file
awk '$3=="exon"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 | cut -f1 -d ';' | \
sed 's/Parent=//g' | sed 's/_T00.//g' | sed 's/chr//g' | \
awk '{print $1,$4-1,$5,".",".",$7,$9}' OFS="\t" | \
sort -b -k1,1 -k2,2n -k3,3n   > Zm-B73-REFERENCE-NAM-5.0.1.canon.exon.bed

#Use bedtools to extract the intron feature
bedtools subtract -a Zm-B73-REFERENCE-NAM-5.0.1.canon.gene.bed -b Zm-B73-REFERENCE-NAM-5.0.1.canon.exon.bed > Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.bed

#Calculate the intron length
awk '{print "intron",$3-$2,$7}' OFS="\t" Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.bed > Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.txt

#use the TE annotation file as the input to generate the TE bed file
#Merge the TE feautures that ovelapped with each other
awk '{print $1,$4,$5,$3}' OFS="\t" Zm-B73-REFERENCE-NAM-5.0.TE.gff3 | sed 's/chr//g' | sort -b -k1,1 -k2,2n -k3,3n > Zm-B73-REFERENCE-NAM-5.0.TE.bed
bedtools merge -i Zm-B73-REFERENCE-NAM-5.0.TE.bed > Zm-B73-REFERENCE-NAM-5.0.TE.merged.bed

#Intersect the intron bed file and TE merged bed file
bedtools intersect -wa -wb -a Zm-B73-REFERENCE-NAM-5.0.1.canon.intron.bed -b Zm-B73-REFERENCE-NAM-5.0.TE.merged.bed > Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.bed

#Calculate the TE length intersect with intron
awk '{print "TE",$10-$9,$7}' OFS="\t" Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.bed > Zm-B73-REFERENCE-NAM-5.0.intron_in_TE.txt

#generate bed file of CDS feature from gff3 file
awk '$3=="CDS"' Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3 | cut -f1 -d ';' | \
sed 's/chr//g' | sed 's/ID=//g' | sed 's/_.00.//g' |\
awk '{print $1,$4-1,$5,".",".",$7,$9}' OFS="\t" | awk '$3!="."' | \
sort -b -k1,1 -k2,2n -k3,3n  > Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS.bed
ml BEDTools
bedtools intersect -wo -a Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS.bed -b Zm-B73-REFERENCE-NAM-5.0.TE.merged.bed > Zm-B73-REFERENCE-NAM-5.0.1.canon.CDS_TE_intersect.bed

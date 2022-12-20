#To calculate the length for different element
gff3=/scratch/yz77862/gene_methylation/gene_annotation/canonical/Zm-B73-REFERENCE-NAM-5.0.1.canon.gff3
for element in UTR exon mRNA gene CDS;do
cut -f1 -d ';' ${gff3} | sed 's/chr//g' | sed 's/five_prime_//g' | sed 's/three_prime_//g'| sed 's/ID=//g' | sed 's/Parent=//g' | sed 's/_T00.//g' | sed 's/_P00.//g' | \
awk '$3=="'${element}'"' | awk '{print $9, $5-$4+1}' OFS="\t" > Zm-B73-REFERENCE-NAM-5.0.1.canon.${element}.length
done

#Read this result in R and then use aggregate function to calculate the length per gene

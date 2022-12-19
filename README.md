Table fo contents:

**1. patterns of gene methylation.**

1.1 Define the gene type

The pipeline for methylation levels calculation in different regions made used of cgmaptools (https://cgmaptools.github.io/cgmaptools_documentation/what-is-cgmaptools.html). The format of methylome file is CGmap file with information of sequence context and methylation level of every covered cytosines.

Calculate gene methylation levels

* ```cgmaptools select region -i NAM.CGmap -r NAM_CDS.bed```   #subset the CGmaps only in CDS region
* ```cgmaptools mtr -i NAM_CDS.CGmap -r NAM_gene.bed```#calculate gene methylated level 
* ```awk``` for gene epiallele classification

Generated core gene file

To define a gene to be core gene, we required a gene to present at 26 genomes. For this end, we made use pangene matrix made in Hufford et al. (2021), an excel file provides interlink for pangene and gene. For each pangene ID, it contains the according geneID for 26 genomes, NA indicates missing syntenic gene and multiple gene ID in one cell indicates tandem duplicated genes. Class indicates the pangene status such as core, near core, dispensibale and private genes. In all analysis, we used the core gene set unless specified.

1.2 Picture the methylation environment around genes(meta gene):
* ```cgmaptools bed2fragreg``` and ```awk``` pipeline (see  ```pre_MFG.sh``` ) are used to generated a tab deliminated file with upstream/downstream 3k and 1.5k within genes in 100bp interval for every genes.
* ```cgmaptools mfg -i```#calculate methylated level 

1.3 The distribution of 100 bp upstream TSS site.

```cgmaptools mtr``` is used to calculate the methylated levels for every 100 bp upstream TSS region.

**2. structural and expression patterns in B73 core genes according to methylation categories.**

2.1 The genetic elements inserted into UM/gbM/teM genes

* ```awk``` command is used to calculate CDS/exon/UTR/intron length for every gene and the summary statistis is done by R base(version 4.10)
* ```bedtools intersect``` is used to indentify the TEs overlapping with genes
*  ```bedtools merge``` is to used to collapse the overlapping TE 
* the customed R script is used to summarize the cumulative TE overlapping with genes.

2.2 TEs from different superfamilies overlapping with genes.

*Use python(ASE/3.19.0-intel-2019b-Python-3.7.4) to process TE gff3 files and make a tab delimited file with chromosome number, start position, end position, TE superfamily and TE family. 
*In shell, we then generated the intron bed files using ```awk``` command
*```bedtools intercept -wa -wb -a -b``` to generate the file with gene overlapping intron

2.3 Proportion of genes with UTR
*```awk``` to generated the UTR bed files
*```bedtools intersect -wa -wb -a gene -b UTR``` 
*Use R to do the summary statistics for insertion frequency for genes with diffferent epiallele status

2.4 Proportion of genes with CDS overlapping repeats

*Use ```awk``` to generated the CDS bed file
*```bedtools intersect``` to intersect CDS and TE files
* In R, use ```aggregate``` to calculate the culmulative length of repeats inserted into each gene.
*  ```merge(,,all.x = T)``` function is used to match geneID and according methylation status, NA is replaced with 0 for ```bedtools``` fail to detect the CDS overlapping with repeats. 
* ```table(length>100,methylation)``` function in R is used to calculate the counts the for genes with different epiallele status.

2.5 Expression pattern for genes with different epiallele status

* To call the tpm values for each gene of all 26 maize lines, we use the pipeline developed in Hufford et. al (2021) with slight modification (see methods). 
* To define the constitutive expressed gene, we require gene to have tpm value above 1 in all investigated 10 tissues (leaf tip, leaf middle, leaf base, root, shoot, ear, anther, tassel, endosperm and embryo), silenced gene to have tpm value all below 1 in 10 tissues and tissue specific to have tpm above in at least one tissue but not all ten. To achieve that, we combined all expression data into one files with gene ID as the identifier. Use ```RowSum(tpm>1)``` from matrixStats package in R to calculate how many tissues have expression for each gene. 
 
**3. epiallele stability.**

3.1 Create the methylation matrix as pangene unit. 

* Use R to create a matrix with genes' cumulative CDS lengths by ```merge(core_pan_ID, CDS_length,by.x=T)``` with for loopping 26 maize genomes. Then use ```rowMedians()``` to calculate the singletons' cumulative CDS lengths for each pangene.
* For the tandem duplicates, we use ```aggregate``` the epialleles together by ```paste0``` function. This will act as a constraint to select genes with similar cumulative CDS lenghts.
* With geneID with similar cumulative CDS, we use ```merge``` to match panID to each geneID to create epiallele matrix. 
* With geneID with similar cumulative CDS, we use ```merge``` to match panID to each geneID to create gene copy matrix. 


3.2 Calculate 
* Divide the epiallele matrix into singleton with all genes to be UM, gbM and teM set.
* NULL all singletons, combine the tandem duplicates as a single list for each set. 
* ```table``` is used to count all UM, gbM and teM is each duplicate list.
* Repeat this process with genes equal and greater than 4.

**4. epialleles switches and gene expression.**

* With geneID with similar cumulative CDS, we use ```merge``` to match panID to each geneID to create 10 TPM matrix for 10 different matrix.
* Subset each TPM matrices for stable UM-gbM, UM-teM and gbM-teM, and make the indication matrix for each epiallele.
* Matrix operation is used to calculate the average TPM values for UM, gbM and teM pangenes in each single unstable set.

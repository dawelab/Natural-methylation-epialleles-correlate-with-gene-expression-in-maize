'''#This script will extract canonical isoforms from a GFF3 gene annotation file'''


#USER INPUT; put location of input and output files here
inPath = '/Users/gent/Desktop/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3'
outPath = '/Users/gent/Desktop/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canon.gff3'


#Open files

inFile = open(inPath, 'r')
outFile = open(outPath, 'w')


geneCount = 0
canonicalCount = 0

		

#Read through entire input file line by line and copy all relevant lines except non canonical transcript isoforms.

for line in inFile:
	try:		
		if line.split('\t')[2] == 'chromosome':
			outFile.write(line) 
		elif line.split('\t')[2] == 'gene':
			outFile.write(line)
			geneCount += 1 #count genes 
		elif line.split('\t')[2] == 'mRNA':
			if 'canonical' in line:
				canonical = True
				canonicalCount += 1 #count canonical isoforms
				outFile.write(line)
			else:
				canonical = False
		else:
			if canonical:
				outFile.write(line) #write all lines related to canonical isoforms
		
	except IndexError: #for header lines, which are not tab delimited
		outFile.write(line)

if geneCount > canonicalCount:
	print('WARNING! MORE GENES THAN CANONICAL ISOFORMS')
	
if canonicalCount < geneCount:
	print('WARNING! LESS GENES THAN CANONICAL ISOFORMS')

	
inFile.close()	
outFile.close()

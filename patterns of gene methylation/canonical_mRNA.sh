
for GENOME in *gff3;do
OUT=${GENOME}_canon.py 
I=${GENOME}
  echo "inPath = '${I}'" >> ${OUT} 
  echo "outPath = '${GENOME}.canon'" >> ${OUT} 
  echo " " >> ${OUT} 
  echo "inFile = open(inPath, 'r')" >> ${OUT} 
  echo "outFile = open(outPath, 'w')" >> ${OUT} 
  echo "geneCount = 0" >> ${OUT} 
  echo "canonicalCount = 0" >> ${OUT} 
  echo " " >> ${OUT} 
  echo "for line in inFile:" >> ${OUT} 
  echo "    try:" >> ${OUT} 
  echo "        if line.split('\t')[2] == 'chromosome':" >> ${OUT} 
  echo "            outFile.write(line)" >> ${OUT} 
  echo "        elif line.split('\t')[2] == 'gene':" >> ${OUT} 
  echo "            outFile.write(line)" >> ${OUT} 
  echo "            geneCount += 1 #count genes" >> ${OUT} 
  echo "        elif line.split('\t')[2] == 'mRNA':" >> ${OUT} 
  echo "            if 'canonical' in line:" >> ${OUT} 
  echo "                canonical = True" >> ${OUT} 
  echo "                canonicalCount += 1 #count canonical isoforms" >> ${OUT} 
  echo "                outFile.write(line)" >> ${OUT} 
  echo "            else:" >> ${OUT} 
  echo "                canonical = False" >> ${OUT} 
  echo "        else:" >> ${OUT} 
  echo "            if canonical:" >> ${OUT} 
  echo "                outFile.write(line) #write all lines related to canonical isoforms" >> ${OUT} 
  echo "    except IndexError: #for header lines, which are not tab delimited" >> ${OUT} 
  echo "        outFile.write(line)" >> ${OUT} 
  echo "if geneCount > canonicalCount:" >> ${OUT} 
  echo "    print('WARNING! MORE GENES THAN CANONICAL ISOFORMS')" >> ${OUT} 
  echo "if canonicalCount < geneCount:" >> ${OUT} 
  echo "    print('WARNING! LESS GENES THAN CANONICAL ISOFORMS')" >> ${OUT} 
  echo "inFile.close()" >> ${OUT} 
  echo "outFile.close()" >> ${OUT} 

done 

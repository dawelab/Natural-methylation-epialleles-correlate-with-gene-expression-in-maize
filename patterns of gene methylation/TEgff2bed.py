#!/usr/bin/env python3
"""

    TEgff2bed - TEgff2bed.py
  
    Copyright: Yibing Zeng
    Contact: Yibing.Zeng@uga.edu

"""

import pandas as pd
from optparse import OptionParser
from pandas import DataFrame

def file(f,filename):
    #Read in the file     
    df = pd.read_csv(f, sep = "\t", skiprows = range(6),header = None, compression = 'gzip')
    
    # names the dataframe
    df.columns = ["chr","method","superfamily","score","start","end","strand","phase","attributes"]
    
    #Extract the TE family information
    family = DataFrame(df.iloc[:,8].str.split("Name=", n =2, expand = True).iloc[:,1].str.split(";", n =4, expand = True)).iloc[:,0]
    
    df['family'] = family
    
    data = df.iloc[:,[0,3,4,2,9]]
    
    data = pd.DataFrame(data)
    
    return(data.to_csv(filename,sep ="\t",index=None,na_rep='NAN', header = None))
    
 
 # ===========================================

    

def main():
    usage = "Usage: .\TEgff2bed.py [-i <TE.gff3.gz>] [-f <filename/path>] \n" \
            "Description: Extract the TE's superfamily and family information and generate bed file \n" \
            "Contact:     Yibing Zeng; YibingZeng@uga.edu\n" \
            "Last Update: 11/04/2021 \n" \
            "Output Ex:\n" \
            "   1      19530    20033    Gypsy_LTR_retrotransposon    xilon_diguus_AC185317_1058 \n" \
            "   1      19663    20114    Copia_LTR_retrotransposon    chr6_D_62813001 \n" \
            "   1      20021    25373    CACTA_TIR_transposon         DTC_ZM00085_consensus \n" \
            "   1      25258    25512    hAT_TIR_transposon           TA_ZM00085_consensus" \
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="input",
                      help= "File name end with gff.gz. "
                            "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-f",  dest="filename",
                      help= "The filename/path given to the output file")
  #
    #
    (options, args) = parser.parse_args()
    #
    file(options.input,options.filename)
#
# ===========================================
if __name__ == "__main__":
    main()
#

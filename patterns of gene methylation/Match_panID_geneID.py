import os
import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.core.frame import DataFrame
data = pd.read_csv("pan_gene_matrix_v3_cyverse.csv")
pan = data['Pan_gene_ID']
data_core = data[data['class']=="Core Gene"]
df['pan'] = data_core['Pan_gene_ID']
pan = data_core['Pan_gene_ID']
def core(NAM,name):
    df_pan = df.iloc[:,NAM]
    new = df_pan.str.split(";", n = 200, expand = True)
    df2 = DataFrame(new)
    df2['pan'] = pan
    df2['copy'] = df2.shape[1] - df2.isna().sum(axis =1) -1
    df3 = DataFrame(df2)
    return(df3.to_csv(name,na_rep='NA',index = None))
list = df.columns.values.tolist()
for i in range(0,26):
    core(i,list[i] + ".core")

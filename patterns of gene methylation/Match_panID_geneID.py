import os
import pandas as pd
os.chdir("/Users/x/Desktop/Data/core_gene")
os.getcwd()
from pandas.core.frame import DataFrame
data = pd.read_csv("/Users/x/Desktop/Data/pan_gene/pan_gene_matrix_v3_cyverse.csv")
pan = data['Pan_gene_ID']
pangeneclass = data['class']
list = data.columns.values.tolist()
def core(NAM,name):
    df_pan = data.iloc[:,NAM]
    new = df_pan.str.split(";", n = 200, expand = True)
    df = DataFrame(new)
    df['pan'] = pan
    df['pangeneclass'] = pangeneclass
    df = DataFrame(df)
    return(df.to_csv(name,na_rep='NA',index = None))
list = data.columns.values.tolist()
for i in range(3,29):
    core(i,list[i] + ".class.txt")

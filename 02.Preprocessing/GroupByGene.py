import pandas as pd
import numpy as np
import re
import sys

inFile = sys.argv[1]
s_table = sys.argv[2]
s_num = re.match(".*/(.*).TNscope",inFile).group(1)

Samples = {}
with open(s_table) as f:
    for line in f:
        s, v = line.strip().split('\t')
        Samples[v] = s

df = pd.read_csv(inFile,sep='\t')

n_df = df.loc[:,["Gene.refGene","GEN[*].AF"]]
n_df.columns = ["Gene","AF"]
n_df.insert(loc=0, column='Sample', value=Samples[s_num])

gbs = n_df.groupby(["Sample","Gene"])
res = gbs['AF'].agg([np.max]).reset_index()
res.columns = ['Sample','Gene','AF']
res = res.sort_values(by=['Gene'], ascending=[True]).reset_index(drop=True)

res.to_csv(f'Sel/{Samples[s_num]}.txt',index=False,sep='\t')

print(f'{s_num}\t{Samples[s_num]}\t{res.shape[0]}')

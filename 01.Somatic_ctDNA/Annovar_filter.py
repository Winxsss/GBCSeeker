#!/usr/bin python
import pandas as pd
import numpy as np
import os
import sys
import re

infile = sys.argv[1]
sample = os.path.splitext(infile)[0]

noUMI = []

df = pd.read_csv(infile,sep='\t')

refGene_requried = df['Func.refGene'].map(lambda x: x in ['exonic','splicing'])
ExonicFunc_required = df['ExonicFunc.ensGene'].map(lambda x : x in ['frameshift_deletion','frameshift-insertion','nonsynonymous_SNV','stopgain'])
bigShift_filter = df['REF'].map(lambda x : len(x) < 10) & df['ALT'].map(lambda x : len(x) < 10)
Exac_required = df['ExAC_ALL'].map(lambda x : x <= 0.001)

def AD_sel(ads):
    s = sample.split('.')[0]
    threshold = 2
    if s in noUMI:
        threshold = 2
    ad = float(ads.split(',')[-1])
    if ad > threshold :
        return True
    else :
        return False

def DP_sel(ads):
    dp = sum(map(lambda x: float(x), ads.split(',')))
    threshold = 10
    if dp > threshold :
        return True
    else :
        return False

AD_required = df["GEN[*].AD"].map(AD_sel)
DP_required = df["GEN[*].AD"].map(DP_sel)

Res = df[refGene_requried & ExonicFunc_required & bigShift_filter & Exac_required & AD_required & DP_required]

def AF_filter(x):
    threshold1 = 0.001
    threshold2 = 0.3
    try:
        if threshold1 < x <= threshold2 :
            return True
        else :
            return False
    except Exception as e :
        return any(map(lambda i : threshold1 < float(i) <= threshold2, re.split(',',x)))

Res.loc[Res["GEN[*].AF"].map(AF_filter)].to_excel(f'{sample}.AD2.xlsx',index=False)
Res.loc[Res["GEN[*].AF"].map(AF_filter)].to_csv(f'{sample}.AD2.txt',index=False,sep='\t')

print(f'## {sample} DONE!')


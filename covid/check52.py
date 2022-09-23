#encoding: utf-8
import os
import pandas as pd
import numpy as np
import sys 
reload(sys) 
sys.setdefaultencoding('utf-8')

os.chdir(r'.')
import tkinter as tk
from tkinter import filedialog as fd
root = tk.Tk()
root.withdraw()
path1 = fd.askopenfilename()
print(path1)

root = tk.Tk()
root.withdraw()
path2 = fd.askopenfilename()
print(path2)

data = pd.read_csv("lims.txt", sep="\t")
dy = pd.read_csv("yundd.txt", sep="\t") 


df = data
dxg = pd.DataFrame(df.groupby([xm, sjdw])[ybtm, sfzh].nunique()).reset_index()
#dxg = pd.DataFrame(df.groupby(xm)[ybtm, sfzh].nunique()).reset_index()
jm_dxg = [1]*dxg[xm].count()
dxg['joinme'] = jm_dxg


dy = dy.iloc[2:,[0,13,1,10,12]]
dy.columns = ['no',xm,'dw','gs','rs']
dxg.columns = [xm,'dw','gs','rs','joinme']
dxg[xm] = dxg[xm].astype(str)
dgy = dxg.merge(dy, on=[xm], how='left')


jm_df=[]
xiang  = {}
for index, dfgrow in df.iterrows():
    xing = dfgrow[xm]
    if xing in xiang:
        jm_df.append(0)
    else:
        xiang[xing] = 1
        jm_df.append(1)

df['joinme'] = jm_df


df[xm] = df[xm].astype(str)
df2 = df.merge(dgy, on=[xm,'joinme'], how='left')
df2 = df2.drop(columns=[u'joinme'])
df2 = df2.replace(np.nan, '', regex=True)
df2.to_csv('hh3.txt', sep='\t')
df2.to_excel('hh3.xlsx')

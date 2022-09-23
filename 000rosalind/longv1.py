!pip install Bio
#!jupyter notebook --NotebookApp.iopub_data_rate_limit=1.0e10

from Bio import SeqIO
import pandas as pd
from io import StringIO
#from pprint import pprint

#The dataset is guaranteed to satisfy the following condition:
#there exists a unique way to reconstruct the entire chromosome
#from these reads by gluing together pairs of reads that 
#overlap by more than half their length.

def lap(fq):
  n = len(fq)
  m = 0.5
  hh = [(0, 0, 0, "npq")] * n
  for i in range(0, n):
    for j in range(i+1, n):
      o = max(x for x in range(len(fq[i])+1) if fq[i].endswith(fq[j][:x]))
      #if o > int(len(fq[i])/2+1):
      if o > len(fq[i])*m:
        hh[i] = i, j, o, "apres"
      o2 = max(x for x in range(len(fq[j])+1) if fq[j].endswith(fq[i][:x]))
      if o2 > len(fq[i])*m:
        hh[j] = j, i, o2, "avant"
  return hh


#tmp walk through the jo cat matrix, i+j-o
def jolympique(df, fqs):
  #len(df.shape[0])
  for i in range(0, len(df.index)):
    j = df.iloc[i][2]
    if j == 0: 
      h = i+1
  j = df.iloc[h][2]
  i = df.iloc[h][1]
  long = fqs[h] + fqs[i][j:]
  while j != 0: 
    p = df.iloc[i][1]
    o = df.iloc[i][2]
    long += fqs[p][o:]
    j = df.iloc[p][2]
    if j == 0:
      print(long,"\n")
      break
#    print(p)
    i = df.iloc[p][1]
    o = df.iloc[p][2]
    long += fqs[i][o:]
#    print(i)
    j = df.iloc[i][2]
  print(long,"\n")


"""    p = df.iloc[i][2]
    if p == "apres":
      long = long+fqs[j][o:]
    elif p == "avant":
      long = fqs[j][:o]+long"""



#if __name__ == "__main__":
fqs = """
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC
>Rosalind_56
ATTAG
ACCTG
"""
f = StringIO(fqs)
arr = []
#arr = fqs.split('\n')
#reads = [i for i in arr if ">" not in i]
for record in SeqIO.parse(f, "fasta"):
  #print(record.id, record.seq)
  arr.append(str(record.seq))
hh = lap(arr)
x = pd.DataFrame(hh)
#print(x)
jolympique(x, arr)

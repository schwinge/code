buffer
>>> from Bio import SeqIO
>>> my_fasta_file = "fuck"
>>> my_sequence = "G"     //?
>>> for record in SeqIO.parse(my_fasta_file, "fasta"):
...     total += record.seq.count(my_sequence)
... 
>>> print total
max0 format

>>> file = open ("haha", "rU")
>>> s = []
>>> for record in SeqIO.parse(file,"fastq"):
...     ph = record.letter_annotations["phred_quality"]
...     s = s+ph
c,n = 0,0
>>> c,n = 0,0
>>> for i,e in enumerate(s):
...     c = c+e
...     print i,e,c/(i+1)
#for i in range(0,len(s)+1):
>>> type(s)
<type 'list'>
>>> type(s[0])
<type 'int'>
>>> file = open ("haha", "rU")
>>> s = []
>>> l,n = 0,0
>>> for record in SeqIO.parse(file,"fastq"):
...     ph = record.letter_annotations["phred_quality"]
...     s = s+ph
...     l = l+1
>>> for i in range(0,len(ph)):
...     c =0
...     for j in range(0,l):
...         c = c+s[j*len(ph)+i]
...     if c < foo*l:
...         n = n+1

>>> for record in SeqIO.parse(file,"fastq"):
...     print record[:1]
... 
ID: Rosalind_0049
Name: Rosalind_0049
Description: Rosalind_0049
Number of features: 0
Per letter annotation for: phred_quality
Seq('G', SingleLetterAlphabet())

trimseq.append(rec[:i-2])
print "Saved %i reads" % len(trimseq)
output_handle=open("C:\work in SF\qtrim_seqs.fastq", "w")
SeqIO.write(trimseq, output_handle, "fastq")
output_handle.close()    #https://www.biostars.org/p/100281/


    :::python
    #!/usr/bin/python2

    import sys
    from Bio import SeqIO

    file = open (sys.argv[1], "rU")
    x = int(sys.argv[2])
    trimseq=[]
    for record in SeqIO.parse(file,"fastq"):
        ph = record.letter_annotations["phred_quality"]
        l = len(ph)
        c = []
        for i,e in enumerate(ph):
            if e < x:
                c.append(0)
            else:
                c.append(1)
        for i in range(0,len(c)-1):
            if int(c[i]+c[i+1])>1:
                m = i
                break
        for i in range(0,len(c)):
            if int(c[len(c)-i-2]+c[len(c)-i-1])==1:
                n = len(c)-i-1
                break
        print m,n
        trimseq.append(record[m:n])
    output_handle=open("check", "w")
    SeqIO.write(trimseq, output_handle, "fastq")  #kinda failed = =
    output_handle.close()









ruby 2.0.0p648 (2015-12-16) [x86_64-linux]
>>> import os
>>> cwd = os.getcwd()
>>> cwd



@click.command()
https://click.palletsprojects.com/en/7.x/


==============================================
ORF Open Reading Frames
#!/usr/bin/python

import sys
#import re

mydict={
    'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
    'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
    'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
    'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
    'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
    'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
    'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'TAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'TAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
    'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'TGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}


def findc(f):
  i = []
  while f:
    c=f.read(3)
    i.append(c)
    if not c:
      break
  return i


def findfv():
  r = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
  f = open(sys.argv[1])
  v=""
  while f:
    c = f.read(1)
    if c in r:
      v=v+r[c]
    if not c:
      break
  return v[::-1]


def finddict(c):
  s=""
  for i in range(0,len(c),1):
    if c[i] in mydict:
      s=s+mydict[c[i]]
  return s
#  s = s.replace("Stop","stop\n")


def shadow(s):
  k=[]
  r=[]
  for i in range(0,len(s),1):
    for j in range(1,len(s[i]),1):
      if s[i][j]=="M":
         k.append(s[i][j:])
  for i in range(0,len(k),1):
    for j in range(1,len(k[i]),1):
      if k[i][j:j+4]=="Stop":
        r.append(k[i][:j])
        break
  return r


if __name__ == "__main__":
#  ss=""
  ss=[]
  for i in range(0,3):
    f = open(sys.argv[1])
    f.read(i)
    c = findc(f)
    s = finddict(c)
#    ss=ss+s
    ss.append(s)

  v = findfv()
  with open("hh","w")as fv:
    fv.write(v)

  for i in range(0,3):
    fv = open("hh")
    fv.read(i)
    cv = findc(fv)
    sv = finddict(cv)
#    ss=ss+sv
    ss.append(sv)
  sk=shadow(ss)
  print '\n'.join(set(sk))



#  es = set(ss.split("\n"))
#  sf=[e.strip("stop") for e in sk if e.startswith("M") and e.endswith("stop")]
#  print '\n'.join(sf)


/M.*\n
#regex doesnt work for MXXMXstop
============================================




=====================================
Locating Restriction Sites

#!/usr/bin/python 
 
 
import sys 
 
 
def findfv(f): 
  r = {'A':'T', 'T':'A', 'G':'C', 'C':'G'} 
  v = "" 
  while f: 
    c = f.read(1) 
    if c in r: 
      v=v+r[c] 
    if not c: 
      break 
  return v 
 
 
def findji(c,v): 
  s=[]   for i in range(4,13,1):     for j in range(0,len(c)-i+1,1):       if c[j:j+i] == v[j:j+i][::-1]: 
        s.append(str(j+1)+" "+str(i)) 
  return s 
 
 
if __name__ == "__main__": 
  f = open(sys.argv[1]) 
  v = findfv(f) 
  c = open(sys.argv[1]).readline().strip("\n") 
  s = findji(c,v) 
  print '\n'.join(s) 
==========================================





===============================================
RNA Splicing
#!/usr/bin/python

import sys
import re


mydict={}


def pro(c):
  s=""
  for i in range(0,len(c)-2,3):
    if c[i:i+3] in mydict:
      s=s+mydict[c[i:i+3]]
  return s


def exon():
  with open(sys.argv[1]) as f:
    y = next(f)
    for l in f:
      line = str(l.strip('\n'))
      if line in y:
        y = y.replace(line,'')
  return y


if __name__ == "__main__":
  c = exon()
  m = pro(c)
  x = re.findall('M.*Stop',m)
  print ''.join(x).strip('Stop')
=================================================

Enumerating k-mers Lexicographically
#!/usr/bin/python

import sys
import itertools

abc=str(sys.argv[1]).replace(" ","")
n=int(sys.argv[2])
kmer=itertools.product(abc,repeat=n)
#kmer=itertools.permutations(abc,n)
for e in kmer:
  print ''.join(e)

===================================================



==============================================

d-score
#!/usr/bin/python

import sys
#python d.py a.d.test a.m.1 > test
i,j = 0,0
d = [line.split() for line in open(sys.argv[1],"r")]
m = [line.split() for line in open(sys.argv[2],"r")]
for i in range(0,len(d)):
        s,ss = 0,0      #for test= 1
        for j in range(0,len(m)):
                ds = float(d[i][0])
                de = float(d[i][1])
                ms = float(m[j][0])
                me = float(m[j][1])     #could be ms
                if ms >= ds and ms <= de:
                        ss += float(m[j][2])
                        if me <= de:
                                s += float(m[j][2])
                if me >= ds and me <= de:
                        ss += float(m[j][2])
        print(s/ss)
==================================================

Longest Increasing Subsequence
numpy set matrix
sn=s[s!=0]
Combinations are emitted in lexicographic sort order    #https://stackoverflow.com/questions/50651619/is-itertools-combinations-deterministic

#!/usr/bin/python

import sys
import itertools

def combi(t,n):
  s=[]
  ss=[]
  for i in range(2,n,1):
    com=itertools.combinations(t,i)
    for e in com:
      if ordertest(list(e)) and len(s)<len(e):
        s=e
      if ins(list(e)) and len(ss)<len(e):
        ss=e
  return s,ss

def ordertest(A):
    return all(A[i] >= A[i+1] for i in range(len(A)-1))

def ins(A):
    return all(A[i] <= A[i+1] for i in range(len(A)-1))

if __name__ == "__main__":
  nn = int(sys.argv[1])
  t = [line.split() for line in open(sys.argv[2],"r")]
  s,ss=combi(t[0],nn)
  print(' '.join(map(str,s)))
  print(' '.join(map(str,ss)))





=================================================================
'''
Charles version
The example from:
Python / OOP / 8: A Simple Game

Extended list of commands:
"say", "examine", "hit", "question"

Extended list of objects:
"goblin", "orc", "human"

For example:
"say hi\n examine goblin\n hit goblin\n hit goblin\n hit goblin\n examine human\n question human\n hit human\n examine human"
'''

## SUPERCLASS
class GameObject:
  class_name = ""
  desc = ""
  objects = {}

  def __init__(self, name):
    self.name = name
    GameObject.objects[self.class_name] = self

  def get_desc(self):
    return self.class_name + "\n" + self.desc

## GOBLIN CLASS
class Goblin(GameObject):
  def __init__(self, name):
    self.class_name = "goblin"
    self.health = 3
    self._desc = " A foul creature"
    super().__init__(name)

  @property
  def desc(self):
    if self.health >=3:
      return self._desc
    elif self.health == 2:
      health_line = "It has a wound on its knee."
    elif self.health == 1:
      health_line = "Its left arm has been cut off!"
    elif self.health <= 0:
      health_line = "It is dead."
    return self._desc + "\n" + health_line
    
  @desc.setter
  def desc(self, value):
    self._desc = value

## HUMAN CLASS
class Human(GameObject):
  def __init__(self, name):
    self.class_name = "human"
    self._desc = "Just a regular human person like you"
    super().__init__(name)
    
  @property
  def desc(self):
    return self._desc
  
  @desc.setter
  def desc(self, value):
    self._desc = value

## ORC CLASS
class Orc(GameObject):
  def __init__(self, name):
    self.class_name = "orc"
    self.health = 5
    self._desc = " A powerful evil creature"
    super().__init__(name)

  @property
  def desc(self):
    if self.health >= 5:
      return self._desc
    elif self.health == 4:
      health_line = "Not even a scratch."
    elif self.health == 3:
      health_line = "You have broken its shield."
    elif self.health == 2:
      health_line = "It has a wound on its knee."
    elif self.health == 1:
      health_line = "Its left arm has been cut off!"
    elif self.health <= 0:
      health_line = "It is dead."
    return self._desc + "\n" + health_line

  @desc.setter
  def desc(self, value):
    self._desc = value

## Input function to interpret the list of commands
def get_input():
  user_input = ""
  while 1:
    try:
      user_input = user_input + input() + "\n"
    except:
      break
  default_input = "say hi\n examine goblin\n question goblin\n hit goblin\n examine goblin\n hit goblin\n examine goblin\n hit goblin\n examine goblin\n examine orc\n question orc\n hit orc\n examine orc\n examine human\n question human\n hit human\n examine human"
  text = user_input if user_input != "\n" else default_input
  comm_list = text.splitlines()
  i = 0
  for line in comm_list:
    command = line.split()
    if len(command) == 1:
      command.append("nothing")
    verb_word = command[0]
    i = i + 1
    print (i,".","You",command[0],command[1])
    if verb_word in verb_dict:
      verb = verb_dict[verb_word]
    else:
      print("Unknown verb {}". format(verb_word))
      return

    if len(command) >= 2:
      noun_word = command[1]
      print (verb(noun_word))
    else:
      print(verb("nothing"))

## Function to react to 'say'
def say(noun):
  return 'You said "{}"'.format(noun)

## Function to react to 'examine'
def examine(noun):
  if noun in GameObject.objects:
    return GameObject.objects[noun].get_desc()
  else:
    return "There is no {} here.".format(noun)

## Function to react to 'hit'
def hit(noun):
  if noun in GameObject.objects:
    thing = GameObject.objects[noun]
    if type(thing) == Goblin:
      thing.health = thing.health - 1
      if thing.health <= 0:
        msg = "You killed this servant of evil!"
      else: 
        msg = "You hit the {}".format(thing.class_name)
    elif type(thing) == Orc:
      thing.health = thing.health - 1
      if thing.health <= 0:
        msg = "You killed this servant of evil!"
      else: 
        msg = "You hit the {}".format(thing.class_name)
    elif type(thing) == Human:
      thing.desc = "You tried hitting this poor person. Now it's mad."
      msg = "Why would you hit a human?"
      
  else:
    msg ="There is no {} here.".format(noun) 
  return msg

## Function to react to 'question'
def question(noun):
  if noun in GameObject.objects:
    creature = GameObject.objects[noun]
    try:
      if creature.health > 0:
        msg = "\"My name is " + creature.name + ".\""
      else:
        msg = "Why questioning the dead..."
    except:
      msg = "\"My name is " + creature.name + ".\""
  else:
    msg ="There is no {} here.".format(noun)
  return msg

# Initialising the program with 1 object per class
goblin = Goblin("Gobbly")
orc = Orc("Orcley")
human = Human("Charles")

verb_dict = {
  "say": say,
  "examine": examine,
  "hit": hit,
  "question": question,
}

get_input()

====================================================

permutation recovery
https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1713794
#Pydna: a simulation and documentation tool for DNA assembly strategies using python
==============================


>>> X = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
>>> Y = [ 0,   1,   1,   0,   1,   2,   2,   0,   1 ]
>>> yx = zip(Y, X)
>>> yx
[(0, 'a'), (1, 'b'), (1, 'c'), (0, 'd'), (1, 'e'), (2, 'f'), (2, 'g'), (0, 'h'), (1, 'i')]
>>> [x for _,x in sorted(zip(Y,X))]
['a', 'd', 'h', 'b', 'c', 'e', 'i', 'f', 'g']
#https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
























#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ggpubr)
library(ggplot2)
df<-read.table(args[1],head=T)
df <- data.frame(df)
ggscatter(df, x = "a", y = "b", 
     cor.coef = TRUE, cor.method = "spearman",
     xlab = " expression", ylab = " expression")



for (j in 1:nrow(x)) {
    for (i in 1:nrow(y)) {
    a <- cor.test(t(x[j, ]), t(y[i, ]), method = "spearman", exact = FALSE ) 
    if(!(is.na(a$estimate) && is.na(a$p.value))){
    if (a$estimate > 0.4 && a$estimate < 1 && a$p.value < 0.001)
        print(rownames(x[j,]))
        break}
    }
    }



x<-a2
for (j in 1:nrow(x)) {
    for (i in (j+1):nrow(x)){
    	for(k in 1:ncol(x)){			
    		if(x[j,k]-x[i,k]<0 || is.na(x[j,k]) || is.na(x[i,k])){
    			s[k]=0
    		}else{
    			s[k]=1
    		}
    	}
    	if (sum(s==0)/length(s)>0.2 && sum(s==1)/length(s)>0.2){
        	print(paste0(rownames(x[j,]),"aand",rownames(x[i,])))
        	}
    }
} 


for (k in 1:nrow(pp)) {
	g=pp$V1[k]
	g2=pp$V2[k]
	meta$gg2 = ifelse(as.numeric(exp_haha[g,])>as.numeric(exp_haha[g2,]),1,0)
	ss=Surv(time, event) ~ gg2
	model <- coxph(ss, data = meta )
	p=summary(model)$coefficient[,"Pr(>|z|)"]
	if (!is.na(p) && p < 0.05){
		print(paste0(g,"\t",g2,"\t",p))	
	}
}

for (k in 1:nrow(pp)) {
	g=pp$V1[k]
	g2=pp$V2[k]
	meta$gg2 = ifelse(as.numeric(exp_haha[g,])>as.numeric(exp_haha[g2,]),'1','0')
	s<-meta$gg2
	if (sum(s==0)/length(s)>0.2 && sum(s==1)/length(s)>0.2){
		print("yes")}else{
		print("no")
	}
}


a3<-c()
for (k in 1:nrow(pp)) {
	g=pp$V1[k]
	g2=pp$V2[k]
	g3=paste(g, g2, sep="|")
	a3<-c(a3,g3)
	meta[g3]= ifelse(as.numeric(exp_haha[g,])>as.numeric(exp_haha[g2,]),1,0)
}
write.table(a3,"",row.names=FALSE,col.names=FALSE,sep="\t")


df <- read.table("ya.dm.tag")
df <- data.frame(df)
colnames(df) <- c ("stage","value")

calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}

p<-ggplot(df, aes(stage, value, fill = stage)) + 
    stat_summary(fun.data = calc_stat, geom="boxplot") +
    theme_classic() +
    ggtitle("All") +
    scale_fill_manual(breaks = c("increase", "decrease"), 
                       values=c("red", "blue")) +
    stat_compare_means( method = "t.test",label.y=6)+
    ylab("log2(FC+1)") +
    xlab("D-Score")
p + My_Theme
ggsave("yargg.pdf")

=======================================================================


import java.util.*;
public class LongestIncreasingSubsequence {
  public static void main( String[] args ) {
    String in = "";
    String[] parts = in.split(" ");
    int[] list = new int[parts.length];
    for(int i = 0; i < list.length; ++i) {
      list[i] = Integer.parseInt(parts[i]);
    }
    
    int[] lis = longestIncreasingSubsequence(list);
    for(int i = 0; i < lis.length; ++i) {
      System.out.print(lis[i] + " ");
    }
    System.out.println();
    int[] lds = longestDecreasingSubsequence(list);
    for(int i = 0; i < lds.length; ++i) {
      System.out.print(lds[i] + " ");
    }
  }
  
  public static int[] longestIncreasingSubsequence( int[] list ) {
    int n = list.length;
    int len = 1;
    int end = 1;
    int[] l = new int[n];
    int[] prev = new int[n];
    l[0] = 1;
    for(int i = 1; i < n; ++i) {
      l[i] = 1;
      for(int j = 0; j < i; ++j) {
        if(list[j] < list[i]) {
          if(l[i] < 1 + l[j]) {
            l[i] = 1 + l[j];
            prev[i] = j;
          }
        }
      }
      if(len < l[i]) {
        len = l[i];
        end = i;
      }
    }
    
    int[] out = new int[len];
    int i = end;
    int added = len-1;
    while(added >= 0) {
      out[added--] = list[i];
      i = prev[i];
    }
    return out;
  }
  
  public static int[] longestDecreasingSubsequence( int[] list ) {
    int n = list.length;
    int len = 1;
    int end = 1;
    int[] l = new int[n];
    int[] prev = new int[n];
    l[0] = 1;
    for(int i = 1; i < n; ++i) {
      l[i] = 1;
      for(int j = 0; j < i; ++j) {
        if(list[j] > list[i]) {
          if(l[i] < 1 + l[j]) {
            l[i] = 1 + l[j];
            prev[i] = j;
          }
        }
      }
      if(len < l[i]) {
        len = l[i];
        end = i;
      }
    }
    
    int[] out = new int[len];
    int i = end;
    int added = len-1;
    while(added >= 0) {
      out[added--] = list[i];
      i = prev[i];
    }
    return out;
  }
}




=============================================

import java.io.File; 
import java.util.Scanner;
import java.util.Formatter;

public class Myclass{
    public static void main(String[ ] args) {
        try{
            Formatter f = new Formatter("f.txt");
            f.format("%s %s %s","1","3","5\r\n");
            f.format("%s %s %s","2","4","6");
            f.close();
            
            File x = new File("f.txt");
            Scanner s = new Scanner(x);
            while(s.hasNext()){
                System.out.println(s.next());
            }
            s.close();
        }catch(Exception a){
            System.out.println("the end");
        }
    }
}
=============================================

import java.util.Iterator;
import java.util.LinkedList;

public class Myclass { 
    public static void main(String[ ] args) {
        LinkedList<String> xx =  new LinkedList<String>();
        xx.add("ss");
        xx.add("sss");
        
        Iterator<String> x = xx.iterator();
        while(x.hasNext()){
            System.out.println(x.next());
        }
    }
}
================================================

import java.util.Iterator;
import java.util.LinkedList;

public class Myclass {
    public static void main(String[ ] args) {
        LinkedList<Integer> ncn = new LinkedList<Integer>();
        ncn.add(1);
        ncn.add(2);
        
        Iterator<Integer> f = ncn.iterator();
        Integer v = f.next();
        f.remove();
        int w = f.next();
        boolean x = f.hasNext();
        System.out.println(v);
        System.out.println(w);
        System.out.println(x);
    }
}
==============================================

import java.util.ArrayList;
import java.util.Collections;

public class Myclass {
    public static void main(String[ ] args) {
        ArrayList<Integer> v = new ArrayList<Integer>();
        v.add(10255);
        v.add(52535);
        v.add(45635);
        v.add(1023455);
        v.add(5438);
        v.add(54634);
        
        Collections.sort(v);
        System.out.println(v);
        Collections.shuffle(v);
        System.out.println(v);
    }
}
===========================================================

import java.util.HashSet;

public class Myclass {
    public static void main(String [ ] args) {
        HashSet<String> purtain = new HashSet<String>();
        purtain.add("merde");
        purtain.add("haha");
        System.out.println(purtain);
    }
}
================================================================

import java.util.HashMap;

public class Myclass {
    public static void main(String[ ] args) {
        HashMap<String, xxx> haha = new HashMap<String, xxx>();
        xxx x = new xxx("haha", "wawa");
        haha.put("not", x);
        System.out.println(haha.get("not").xx);
    }
}

public class xxx {
    public String xx = "intiate";
    public String hh = "aoe";
    
    xxx (String xx, String hh) {
        this.xx = xx;
        this.hh =hh;
    }
}
================================================================

import java.util.LinkedList;

class Myclass {
    public static void main(String[ ] args){
        LinkedList<Integer> ll = new LinkedList<Integer>();
        ll.add(100);
        for (Integer i: ll){
            System.out.println(i);
        }
    }
}
=======================================================

import java.util.ArrayList;

public class Myclass {
    public static void main(String[ ] args) {
        ArrayList<String> trois = new ArrayList<String>();
        trois.add("eins");
        trois.add("zwei");
        System.out.println(trois);
        trois.add("drei");
        trois.add("vier");
        trois.add("fuenf");
        trois.add("sechs");
        trois.add("sieben");
        trois.add("acht");
        trois.add("neun");
        trois.add("zehn");
        System.out.println(trois.get(4));
        trois.remove("vier");
        System.out.println(trois.contains("vier"));
        trois.clear(); 
        System.out.println(trois.size());
    }
}
=============================================================

public class MyClass {
    public static void main(String[ ] args) {
        try{
            int k = 1/0;
            k = k;
        } catch (ArithmeticException x){
            System.out.println("0!");
        }
        try{
            System.out.println("test");
            Thread.sleep(1000);
        } catch (InterruptedException e){
            System.out.println("schlafen???");
        }
    }
}
============================================================

class xx implements Runnable {
    public void run() {
        System.out.println("xxx");
    }
}

class Myclass {
    public static void main(String [ ] args) {
        Thread x = new Thread(new xx());
        x.start();
    }
}
-----------------------------
class xx extends Thread {
    public void run() {
        System.out.println("test");
    }
}

class Myclass {
    public static void main(String[ ] args) {
        xx neu = new xx();
        neu.start(); 
    }
}
========================================================

enum Farbe {
    rot, gruen, blau, weiss, schwarz, gelb, braun;
}


        
      
        
public class Myclass {
    public static void main(String[] args) {
        Farbe a = Farbe.rot;
        Farbe b = Farbe.braun;
        Farbe c = Farbe.blau;
        Farbe d = Farbe.gelb;
            switch(b) {
                case rot:
                    System.out.println("nur");
                    break;
                case braun:
                case blau:
                    System.out.println("Noun");
                    break;
                default:
                    System.out.println("!");
                    break;
            }
    }
}
====================================================

class Animal {
    String name;
    Animal(String n) {
        name = n;
    }
}

class MyClass {
    public static void main(String[ ] args) {
        Animal a1 = new Animal("Robby");
        Animal a2 = new Animal("Robby");
        System.out.println(a1 == a2);
        System.out.println(a1 + "  " + a2);
    }
}
=========================================================

class etwas{
    String haha;
    etwas(String lol){
        haha = lol;
        nichts x = new nichts();
        x.denke();
    }
    
    class nichts{
        public void denke(){
            System.out.println(haha+" ich habe keine Achnung");
        }
    }
}

public class Myclass{
    public static void main(String[] args){
        etwas haha = new etwas("wakaka");
    }
}

==========================================================

class Tier {
    public void Farbe() {
        System.out.println("schwarz");
    }
}


public class Myclass {
    public static void main(String[] args) {
        Tier weiss = new Tier() {
            @Override public void Farbe() {
                System.out.println("weiss");
            }
        };
        Tier suess = new Tier();
        weiss.Farbe();   
        suess.Farbe();   
    }
} 
========================================================

class Machine {
    public void start() {
        System.out.println("Starting...");
    }
}

class Program {
    public static void main(String[ ] args) {
        Machine m = new Machine() {
            @Override public void start() {
                System.out.println("Wooooo");
            }
        };
        m.start();
    }
}
==============================================================

class Natur {
    public void Farbe() {
        System.out.println("schwarz");
    }
}

class Baum extends Natur {
    public void Farbe() {
        System.out.println("gruen");
    }
}

class  Blaum extends Natur {
    public void Farbe() {
        System.out.println("rot");
    }
}

public class Myclass {
    public static void main(String[] args) {
        Natur neu = new Baum();
        Natur neuneu = new Blaum();
        neu.Farbe();
        neuneu.Farbe();
    }
}
================================================================

class haha{
    haha(){
        sayhaha();
    }
    
    private void sayhaha(){
        System.out.println("Private haha");
    }
/*    public void sayhaha() {
        System.out.println("Public haha");
    }*/
}

public class privatemethodsinheritance extends haha{
    public static void main(String[] args){
        privatemethodsinheritance p = new privatemethodsinheritance();
    
        System.out.println("main :-D");
    }
}
==================================================================

class A {
    int num1 = 50; 
    public A() {
        System.out.println("New A");
    }
}
class B extends A {
    int num1 = 100; 
    public B() {
        System.out.println("New B");
        System.out.println(num1); 
        System.out.println(super.num1); 
    }
}

class Program {
    public static void main(String[ ] args) {
        B obj = new B();
    }
}
=================================================================

class A {
    public A() {
        System.out.println("I'm A");
    }
}

class B extends A {
    public B() {
        System.out.println("I'm B");
    }
}

public class grant {
    public static void main(String[] args) {
        B b = new B();
    }
}
========================================================

class Tier {
    protected int Fuss;
    public void essen() {
        System.out.println("kuerchen");
    }
}

class Hund extends Tier {
    Hund() {
        Fuss = 4;
    }
}

public class Myclass {
    public static void main(String[] args) {
        Hund meineHund = new Hund();
        meineHund.essen();
    }
}
=============================================================

public class bored {
    public static void gelangweilt() {
        System.out.println("shut ur maul!");
    }
}

public class Myclass {
    public static void main(String[] args){
        bored.gelangweilt();
    }
}
=============================================================

public class power {
    public static double POW = 2;
    power(){
        POW = Math.pow(POW,2);
    }
}

public class Myclass {
    public static void main(String[] args){
        power p1 = new power();
        power p2 = new power();
        power p3 = new power();
        System.out.println(power.POW);
    }
}
=======================================================

public class ceilingfloor {
    public static void main(String args[]) {
        double d = 5.02562;
        double c = Math.ceil(d);
        double f = Math.floor(d);
        System.out.println(c+" "+f);
    }
}
====================================================

public class Myclass {
    public static void main(String args[]) {
        Auto m,s;
        m = new Auto("gruen");
        neustreichen(m);        
        s = new Auto("blau");
        s.setFarbe("roet");
        System.out.println(m.getFarbe()+s.getFarbe());
    }
    
    static void neustreichen(Auto x) {
        x.setFarbe("weiss");
    }
}

public class Auto {
    private String Farbe;

    Auto(String n) {
        this.Farbe = n;
    }
    
    public String getFarbe() {
        return Farbe;
    }
    
    public void setFarbe(String n) {
        this.Farbe = n;
    }
}
===============================================================

public class Myclass{
    public static void main(String args[]){
        int x = 5;
        adtwo(x);
        System.out.println(x);
    }
    static void adtwo(int x){
    x = x+2;
    }
}
====================================================

public class tier{
    int age;
    double weight;
    String feet;
    
    void beep(){
        System.out.println("wuwuwu");
    }
}

class myclass{
    public static void main(String[] args){
        tier cabaka = new tier();
        tier insekt = new tier();
        insekt.feet = "multiple";
        cabaka.beep();
        
    }
}
==============================================

public class Tier{
    void katze(){
        System.out.println("meow moew moew-mo");
    }
}

class myclass{
    public static void main(String[] args){
        Tier susskatz = new Tier();
        susskatz.katze();
    }
}
=============================================

class sum{
    static int sum(int x, int y){
        return x + y;
    }
    public static void main(String[] args){
        int x = sum(2,5);
        System.out.println(x);
    }
}
==============================================

class hi{
    static void hi(String name){
        System.out.println("hi " + name);
    }
    public static void main(String[] args){
        hi("haha");
        hi("double");
    }
}
=============================================






buffer
>>> from Bio import SeqIO
>>> my_fasta_file = "haha"
>>> my_sequence = "G"     
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
...     if c < haha*l:
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
output_handle.close()    #https:


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




[anqi@node10 fastx_bin]$ ruby -v
ruby 2.0.0p648 (2015-12-16) [x86_64-linux]
>>> import os
>>> cwd = os.getcwd()
>>> cwd
'/lustre/anqi/workzone/tools/fastx_bin'



@click.command()
https:


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

##longtime havent used class(for definiting the elves~)




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
Combinations are emitted in lexicographic sort order    #https:

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
#too long, might use java

#done need shebang #!
[anqi@node11 rn32]$ javac haha.java    #generate HelloWorld.class 
[anqi@node11 rn32]$ java HelloWorld
Hello, World
------------




====================================================


for i in $(ls ../*bed); do
    e=${i%%.bed}
    ii=${e##*/}
    bedtools intersect -a ${i} -b fpkm572 -wb \
    | awk '{print$7,log($9/($8+0.01))/log(2)}'|sort|uniq| \
    sed $'s/\ /\t/g'  > ${ii}.fpkm572;
done


Tool:
https://github.com/ay-lab/FitHiChIP/blob/master/Imp_Scripts/PeakInferHiChIP.sh 


Methodology:

Peaks were called using MACS28 (-q 0.01 --nomodel --extsize 147) //https://media.nature.com/original/nature-assets/nmeth/journal/v15/n3/extref/nmeth.4583-S1.pdf
MACS2 was then run on the BED file using the no model and extsize 147 parameters and an FDR cutoff of 1%.  //paper: https://www.nature.com/articles/nmeth.3999
macs2 callpeak -t k27.vp -n descvp -q 0.01 --nomodel --extsize 147 -g mm //default hs
Dangling-end and self-ligation reads from the GM12878 Smc1a HiChIP HiC-Pro output were combined and processed to be compatible as a MACS2 BED file input.  //SCparis& DEpairs
https://www.biostars.org/p/300640/
Inferring peaks from HiChIP data: For the purposes of calling 1D peaks from HiChIP
data when ChIP-seq peaks are not avialable, we have tested different combinations of following four sets of reads generated from the HiC-pro pipeline: 1) dangling end (DE), 2) self-cycle (SC), 3) re-ligation (RE), and 4) CIS ‘short-range’ (< 1 Kb) valid (V) reads (after duplicate removal). //https://www.biorxiv.org/content/biorxiv/early/2018/09/10/412833.full.pdf
PLAC-seq (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5143423/)
In situ Hi-C or PLAC-seq contact maps were visualized using Juicebox19 370 after removing all trans reads and cis reads pairs span less than 10kb. //https://www.biorxiv.org/content/biorxiv/early/2016/09/09/074294.full.pdf



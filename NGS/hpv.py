test2pl.sh
reference_fasta=hg19.fa
fq=$1
samplename=$2
#rgstr='@RG\tID:\tPL:\tSM:'
ngm=ngm/
tvc=tvc/




time gatk-launch HaplotypeCaller \
 -R $reference_fasta \
 -I  $ngm$2.bam\
 -O  $tvc$2gatk.g.vcf && echo "** gvcf done **"
#joint calling
time gatk-launch GenotypeGVCFs \
 -R $reference_fasta \
 -V $tvc$2gatk.g.vcf \
 -O $tvc$2gatk.vcf  && echo "** vcf done **"

bash test2pl.sh x test2.sorted (seems gatk cant set multithreads)}

tvcassembly-->Long Indel Assembly; tvcutils-->unify_vcf)
vcf-->bgzip(zip),tabix(index)
tvc -r xxx -b xxx -o xxx (cd hpvvc/tvc)
tvc -r hg19.fa -b test.sorted.bam -o test.vcf
echo -e '@RG\tID:id\tPL:iontorrent\tLB:library\tSM:sample' > rg
samtools view -H test.sorted.bam | cat - rg > header
find HPV_fastq -name '*026*' -ls

coverage_to_regions.py ref.fa.fai 50000 >ref.fa.50k.regions
freebayes-parallel ref.fa.50k.regions 36 -f ref.fa aln.bam >results.vcf
scripts hg19.fa.fai 50000 > hg19.fa.50k.regions
freebayes-parallel hg19.fa.50k.regions 20 -f hg19.fa test2.sorted.bam > test2freebaysparallel.vcf


paper{
  https://www.sciencedirect.com/science/article/pii/S0014480016302970?via%3Dihub
  https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-264
  https://www.sciencedirect.com/science/article/pii/S1872497316301417?via%3Dihub
  https://www.sciencedirect.com/science/article/pii/S2405852115000051#bib34
  https://www.sciencedirect.com/science/article/pii/S0009898114001521
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5584341/
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3313698/
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5215284/
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753459/
}
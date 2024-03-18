===================HG00419(WGS CRAM I)
kourami-0.9.6me/scripts/alignAndExtract_hs38.sh -d kourami-0.9.6me/db -r kourami-0.9.6me/resources/hs38NoAltDH.fa HG00419 currentworkzone/HG00419.a.sort.bam
export PATH=xx
java -jar kourami-0.9.6me/target/Kourami.jar -d kourami-0.9.6me/db -o HG00419 currentworkzone/HG00419_on_KouramiPanel.bam
=================NA06985 (WES)
biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/bamsort SO=coordinate I=datas/NA06985.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram O=NA06985.e.a.sort.bam index=1 indexfilename=A06985.e.a.sort.bam.bai inputthreads=40 outputthreads=40 fixmates=1 reference=kourami-0.9.6me/resources/hs38DH.fa inputformat=cram outputformat=bam(NAbai!)
WARNING: SAM header designates more than one PG tree root by PP tags.?
currentworkzone/NA06985.e.a.sort.bam
kourami-0.9.6me/scripts/alignAndExtract_hs38.sh -d kourami-0.9.6me/db -r kourami-0.9.6me/resources/hs38NoAltDH.fa NA06985 currentworkzone/NA06985.e.a.sort.bam(bamUtil)
currentworkzone/NA06985_on_KouramiPanel.bam
java -jar kourami-0.9.6me/target/Kourami.jar -d kourami-0.9.6me/db -o NA06985 currentworkzone/NA06985_on_KouramiPanel.bam
slightly difference, try hs38DH
vim ~/.bash_profile
source ~/.bash_profile(or restart)
kourami-0.9.6me/scripts/alignAndExtract_hs38.sh -d kourami-0.9.6me/db -r kourami-0.9.6me/resources/hs38DH.fa NA06985 currentworkzone/NA06985.e.a.sort.bam**the sam bam e.a, cover problem on NA06985_on_KouramiPanel.bam


==== Aligning WGS Reads to the Human Reference Genome:
bwa mem -t 16 -B 4 -O 6 -E 1 -M hs38DH HLA/kourami-0.9.6/resources/hs38DH.fa HLA/kourami-0.9.6/data/AK1/test1m.r1.rq HLA/kourami-0.9.6/data/AK1/test1m.r2.rq | samtools view -1 - > test1m.bam -X: not @RG etc. full of errors
bwa mem -t 16 HLA/kourami-0.9.6/resources/hs38DH.fa HLA/kourami-0.9.6/data/AK1/test1m.r1.rq HLA/kourami-0.9.6/data/AK1/test1m.r2.rq | samtools view -1 - > test16m.bam***[W::bseq_read] the 1st file has fewer sequences. &[gzclose] buffer error
==== HLA Read Extraction and Realignment
HLA/kourami-0.9.6/scripts/alignAndExtract_hs38DH.sh -d HLA/kourami-0.9.6/db -r HLA/kourami-0.9.6/resources/hs38DH.fa test1m HLA/currentworkzone/test16m.bam**missing bamUtil(PATH=tools/bamutil:$PATH**tools/bamutil:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/lib64/mpich-3.2/bin:/usr/java/default/bin:/lustre/apps/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin make install INSTALLDIR=tools/bamutil)**missing Kourami Panel sequences    : HLA/kourami-0.9.6/custom_db/3.33.0/All_FINAL_with_Decoy.fa.gz, so original db 3.24.0
==== Running Kourami
java -jar HLA/kourami-0.9.6/target/Kourami.jar -d HLA/kourami-0.9.6/db -o test1m HLA/currentworkzone/test1m_on_KouramiPanel.bam
#https://github.com/gt1/libmaus2


2. Local realignment around known indels by GATK.
`java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing`
3. Recalibrate base quality scores using known SNPs by GATK.
`java $jvm_args -jar GenomeAnalysisTK.jar -T BaseRecalibrator -nt 1 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R $reference_fasta -o $recal_data.table -I $bam_file -knownSites $known_snps_from_dbSNP142`
`java $jvm_args -jar GenomeAnalysisTK.jar -T PrintReads -l INFO -R $reference_fasta -o $recalibrated_bam -I $bam_file -BQSR $recal_data.table --disable_bam_indexing`







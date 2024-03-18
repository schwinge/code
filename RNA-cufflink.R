library(cummeRbund)
cuff<-readCufflinks()
Creating database /mnt/hwt2_data1/anqi/diff_out/cuffData.db    
csDensity(genes(cuff))
dev.off()
csScatter(genes(cuff), "y", "a")
csVolcano(genes(cuff), "y", "a")  #unlike python #p-adjust maybe
csVolcano(genes(cuff_data_1),'S1','S2', alpha=0.05, showSignificant=T) 
mygene <- getGene(cuff, "Tpp2")
expressionBarplot (mygene)
expressionBarplot(isoforms (mygene))


gene_diff_data <- diffData(genes(cuff))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sig_gene_data)
length(diffGeneIDs)
diffGeneIDs <- getSig(cuff, level="genes",alpha=0.25)  #chosen, Q:"XLOC_000046" "XLOC_000051"
length(diffGeneIDs)
diffGeneIDs <- getSig(cuff, level="genes",alpha=0.1)
length(diffGeneIDs)
diffGeneIDs <- getSig(cuff, level="genes",alpha=0.2)
length(diffGeneIDs)
allgeneFPKMs<-fpkm(genes(cuff))
diffGeneIDs <- getSig(cuff, level="genes",alpha=0.25)
siggeneFPKMs<-allgeneFPKMs[allgeneFPKMs$gene_id %in% diffGeneIDs,]
nrow(siggeneFPKMs)
siggeneFPKMs<-gene_diff_data[gene_diff_data$gene_id %in% diffGeneIDs,] #http://seqanswers.com/forums/archive/index.php/t-18357.html
nrow(siggeneFPKMs)
diffgenes<-getGenes(cuff,diffGeneIDs) 	#http://seqanswers.com/forums/showthread.php?t=18357, thomas doktor
names<-featureNames(diffgenes) #
row.names(names)=names$tracking_id
diffGenesNames<-as.matrix(names)
diffGenesNames<-diffGenesNames[,-1]
diffGenesData<-diffData(diffgenes)
row.names(diffGenesData)=diffGenesData$gene_id
diffGenesData<-diffGenesData[,-1]
diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")
head(diffGenesOutput)

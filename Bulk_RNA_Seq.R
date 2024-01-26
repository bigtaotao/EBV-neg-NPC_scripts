
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

TCGA_HNSC <- read.csv("./TCGA_HNSC_data_biopsy_site_rename.csv")
paracancerous_tissue <- read.csv("./GSE134886/GSE134886_paracancerous_tissue_count.csv")
nasopharyngeal_mucosa <- read.csv("./GSE118719/GSE118719_4_nasopharyngeal_mucosa.csv")

mouse_organoid <- read.csv("./mouse_normal_organoid.csv")
mouse_organoid <- mouse_organoid[!duplicated(mouse_organoid$X),]
dim(mouse_organoid)

all_count <- merge(TCGA_HNSC,paracancerous_tissue,by.x="X",by.y= "X")
all_count <- merge(all_count,nasopharyngeal_mucosa[,c(1,3:6)],by.x="X",by.y= "Geneid")
nrow(all_count)

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(DESeq2)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})

load(file="./ebg_GRCh38.RData")
load(file="./txdb_GRCh38.RData")

all_count$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = as.character(all_count$X),
						keytype ="ENSEMBL",
						column ="SYMBOL",
						multiVals="first") 
library("nichenetr")
aa <- as.character(all_count$symbol) %>% convert_human_to_mouse_symbols()
all_count$mouse_symbol <- as.character(aa)
all_count_1 <- all_count[!duplicated(all_count$mouse_symbol),]
nrow(all_count_1)
all_count_1 <- na.omit(all_count_1)
nrow(all_count_1)

countdata <- merge(all_count_1,mouse_organoid,by.x="mouse_symbol",by.y= "X")
rownames(countdata) <- countdata$mouse_symbol
countdata_mix <- countdata[,c(3:20,34:40,42:44)]
countdata_mix <- countdata[,c(3:40,42:44)] 

colnames(countdata_mix)
id <- colnames(countdata_mix)
species <- c(rep('human',25),rep('mouse',3))
database <- c(rep('TCGA',18),rep('GSE134886',3),rep('GSE118719',4),rep('cclab',3))
biopsy_site <- c(rep('hsa_Base_of_tongue',2),rep('hsa_Floor_of_mouth',3),rep('hsa_Larynx',12),rep('hsa_Mouth',1),rep('hsa_Nasopharynx',7),rep('mus_Nasopharynx',3))
site <- c(rep('other',18),rep('nasoph',10))
countdata_mix_df <- data.frame(id,species,database,biopsy_site,site)
dds_mix <- DESeqDataSetFromMatrix(countData = countdata_mix, colData = countdata_mix_df, design = ~ biopsy_site)

DESeq_counts_mix <- DESeq(dds_mix)
normalized_DEseq_mix <- counts(DESeq_counts_mix, normalized=TRUE)

countdata <- merge(all_count_1,mouse_tissue,by.x="mouse_symbol",by.y= "X")
rownames(countdata) <- countdata$mouse_symbol
countdata_mix <- countdata[,c(3:20,34:40,42:43)]
countdata_mix <- countdata[,c(3:40,42:43)]

colnames(countdata_mix)
id <- colnames(countdata_mix)
species <- c(rep('human',25),rep('mouse',2))
database <- c(rep('TCGA',18),rep('GSE134886',3),rep('GSE118719',4),rep('cclab',2))
biopsy_site <- c(rep('hsa_other',18),rep('hsa_nasoph',7),rep('mus_nasoph',2))
site <- c(rep('other',18),rep('nasoph',9))
countdata_mix_df <- data.frame(id,species,database,biopsy_site,site)
dds_mix <- DESeqDataSetFromMatrix(countData = countdata_mix, colData = countdata_mix_df, design = ~ biopsy_site)

DESeq_counts_mix <- DESeq(dds_mix)
normalized_DEseq_mix <- counts(DESeq_counts_mix, normalized=TRUE)

#limma remove batch effect
library("ggplot2")
pr <- prcomp(log(normalized_DEseq_mix+1,2))
temp <- data.frame(pr$rotation)
load <- data.frame(component=rownames(temp),PC1=temp$PC1,PC2=temp$PC2)
load1 <-cbind(load,countdata_mix_df)
pca.var <- pr$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1) 
p1 <- ggplot(load1, aes(x = PC1, y = PC2, color=biopsy_site,shape = species)) +
geom_point(size =2) +
geom_text(aes(label=id))+
xlab(paste0("PC1: ", pca.var.per[1], "% variance")) +
ylab(paste0("PC2: ", pca.var.per[2], "% variance")) +
labs(title="before reduce batch effect by limma") +
theme_classic() 

limma_edata <- limma::removeBatchEffect(log(normalized_DEseq_mix+1,2), batch = countdata_mix_df$database)

pr_1 <- prcomp(limma_edata)
temp <- data.frame(pr_1$rotation)
load <- data.frame(component=rownames(temp),PC1=temp$PC1,PC2=temp$PC2)
load2 <-cbind(load,countdata_mix_df)

pdf("mouse_organoid on human nolabel.pdf")
ggplot(load2, aes(x = PC1, y = PC2, color=biopsy_site,shape = species)) +
geom_point(size =5) +
xlab(paste0("PC1: ", pca.var.per[1], "% variance")) +
ylab(paste0("PC2: ", pca.var.per[2], "% variance")) +
labs(title="before reduce database batch effect by limma") +
theme_classic()
dev.off()



********************************************************primary_TMP vs normal******************************
********************************************************primary_TMP vs normal******************************
********************************************************primary_TMP vs normal******************************
********************************************************primary_TMP vs normal******************************

#R3.6

project <- c("P_ctrl_VS_normal_organoid")
deal_design <- c("P_ctrl","normal_organoid")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- "./P_ctrl_vs_nomal_organoid/"
sample_sampletable.path <- "./star_mm10_out"
load(file="./ebg_mm10.RData")
load(file="./txdb_mm10.RData")  


deal_case <- deal_design[1]
dela_control <- deal_design[2]
deal_names <- paste(deal_case,dela_control,sep="_v_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"00_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"01_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"01_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"01_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"02_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"03",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"04",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"05",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"06",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"07",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"08",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"09",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(DESeq2)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})

if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
	load(file="./ebg_mm10.RData")
	load(file="./txdb_mm10.RData")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
	load(file="./ebg_GRCh38.RData")
	load(file="./txdb_GRCh38.RData")
}


indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sample_sampletable_ctrl_normal_organoid.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",   
singleEnd=FALSE,  
ignore.strand=TRUE,
fragments=FALSE) 
setwd(file_path)
save(se,file=se.save_RData)

setwd(file_path)
getwd()

colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)

rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(vsd$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_vsd)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()

exons<-as.data.frame(ebg)
exon_length<-data.frame(gene_name=c(0),length=c(0))
for(i in 1:nrow(assay(se))){
  exonsfromgene<-exons[which(exons$group==i),]
  length<-sum(exonsfromgene$width)
  gene_name<-exonsfromgene[1,2]
  exon_length<-rbind(exon_length,c(gene_name,length))
}
exon_length<-exon_length[-1,]
whole_length<-as.data.frame(round(colSums(assay(dds))/1e6,1))
countM<-as.data.frame(assay(se))
FPKM<-countM
for(i in 1:nrow(FPKM)){
  gene<-rownames(FPKM[i,])
  length<-as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
  for(j in 1:ncol(FPKM)){
    FPKM[i,j]<-FPKM[i,j]/(whole_length[j,]*length*0.001) 
  }
}
fpkmToTpm <- function(FPKM){
    exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- sampleTable$sample
colnames(count) <- sampleTable$sample
count_and_tpm <- cbind(count,tpm) 

#DESeq2 normalize
colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","P_ctrl","normal_organoid"))    #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)
write.csv(tpm,file=tpm_csv)


dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)


res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

downres_1  <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

KEGG_pro_enhan_up <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")
KEGG_pro_enhan_down <- setReadable(KEGGdownres_1, "org.Mm.eg.db", keyType="ENTREZID")

write.csv(KEGG_pro_enhan_up, file = KEGGupres_1_file.csv)
write.csv(KEGG_pro_enhan_down, file = KEGGdownres_1_file.csv)

pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()


pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()

##GSEA

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_h -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c1 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c2 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c3 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c4 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c5 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c6 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c7 -gui false



$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3E $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

all_gene <- read.csv("./Primary_TMP_v_normal_organoid_heatmap_data.csv")
head(all_gene)
all_gene <- all_gene[!duplicated(all_gene$X),]
all_gene <- na.omit(all_gene)
rownames(all_gene) <- all_gene[,1]
all_gene <- all_gene[,-1]

colnames(all_gene)
all_gene <- all_gene[,-c(7:9)]
zscore <- t(apply(all_gene, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
head(zscore)

library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
source("./my_code/MyBestFunction_scRNA.R")
source("./my_code/Pseudo_CNV_series.R")

SeuratObject <- CreateSeuratObject(counts = zscore, project = "WXD")
SeuratObject$orig.ident <-factor(SeuratObject$orig.ident,levels = c("normal","ctrl"))
gene <- rownames(zscore)

pdf("NPC_marker_heatmap.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Bmi1","Tubb2a","Tubb2b","Muc5ac","Krt5","Krt13","Krt14","Krt16","Notch1","Egr2","Cd44","Cav1","Cd109"))
dev.off()


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3G $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3G $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig3G $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


**************************GSE13597
library("GEOquery")
gse <- getGEO("GSE13597",GSEMatrix = TRUE, destdir = ".",getGPL = T, AnnotGPL = T)
exprs <- exprs(gse[[1]])
fdata <- fData(gse[[1]])
explan <- data.frame(exprs)
explan$ID <- fdata$ID
explan$symbol <- fdata$"Gene symbol"
explan$gene_ID <- fdata$"Gene ID" 
express_2 <- explan
express_2$symbol <- data.frame(sapply(express_2$symbol,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$gene_ID <- data.frame(sapply(express_2$gene_ID,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$rowMeans <- apply(express_2[,1:28],1,function(x) mean(as.numeric(x),na.rm = T))

express_2 <- express_2[order(express_2$rowMeans,decreasing = T),]
express_na <- na.omit(express_2)
countdata <- express_na[express_na$symbol !="" ,]
explan_final <- express_na[!duplicated(express_na$symbol),]
colnames(explan_final) <- c("Normal1","Normal2","Normal3","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25","array_ID","Gene_symbol","gene_ID","rowMeans")
rownames(explan_final) <- explan_final[,30]
count_anno <- explan_final[,1:28]
write.csv(count_anno, "GSE13597_count_anno.csv")
count_anno_without_arrayID <- explan_final[,1:31]
write.csv(count_anno_without_arrayID, "GSE13597_count_anno_without_arrayID.csv")

group_list <- factor(c(rep('Normal',3),rep('NPC_tumor',25))) 
table(group_list)

library(limma)
library(data.table)
rt <- count_anno

rt <- as.matrix(rt)

dimnames=list(rownames(rt),colnames(rt))

data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data=avereps(data)
data<-as.data.frame(data)
nrow(data)
data=data[rowMeans(data)>0.5,]

design=model.matrix(~factor(group_list))
design
fit=lmFit(data,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)

all_DEG <- merge(count_anno_without_arrayID,DEG, by = "row.names", all=FALSE)
write.csv(all_DEG, "GSE13597_all_DEG.csv")



library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")
library("dplyr")
library("biomaRt")
library("remotes")
library("nichenetr")

all_TF <- all_DEG
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
write.csv(matrix,"GSE13597_all_DEG.csv")



********************************************GSE12452

library("GEOquery")
gse <- getGEO("GSE12452",GSEMatrix = TRUE, destdir = ".",getGPL = T, AnnotGPL = T)
exprs <- exprs(gse[[1]])
fdata <- fData(gse[[1]])
explan <- data.frame(exprs)
explan$ID <- fdata$ID
explan$symbol <- fdata$"Gene symbol"
explan$gene_ID <- fdata$"Gene ID" 
head(explan)
express_2 <- explan
express_2$symbol <- data.frame(sapply(express_2$symbol,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$gene_ID <- data.frame(sapply(express_2$gene_ID,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$rowMeans <- apply(express_2[,1:41],1,function(x) mean(as.numeric(x),na.rm = T))
head(express_2)
express_2 <- express_2[order(express_2$rowMeans,decreasing = T),]
express_na <- na.omit(express_2)

explan_final <- express_na[express_na$symbol !="" ,]

explan_final <- express_na[!duplicated(express_na$symbol),]
colnames(explan_final) <- c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Normal7","Normal8","Normal9","Normal10","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25","NPC26","NPC27","NPC28","NPC29","NPC30","NPC31","array_ID","Gene_symbol","gene_ID","rowMeans")
rownames(explan_final) <- explan_final[,43]
count_anno <- explan_final[,1:41]
write.csv(count_anno, "GSE12452_count_anno.csv")
count_anno_without_arrayID <- explan_final[,1:44]
write.csv(count_anno_without_arrayID, "GSE12452_count_anno_without_arrayID.csv")

group_list <- factor(c(rep('Normal',10),rep('NPC_tumor',31)), levels=c("Normal","NPC_tumor")) 
table(group_list)

library(limma)
library(data.table)
rt <- count_anno
rt <- as.matrix(rt)

dimnames <- list(rownames(rt),colnames(rt))
data <- matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data <- avereps(data)
data <- as.data.frame(data)
nrow(data)
data <- data[rowMeans(data)>0.5,]

design <- model.matrix(~factor(group_list))
design
fit <- lmFit(data,design)
fit <- eBayes(fit)
DEG <- topTable(fit,coef=2,n=Inf)
head(DEG)

all_DEG <- merge(count_anno_without_arrayID,DEG, by = "row.names", all=FALSE)
write.csv(all_DEG, "GSE12452_all_DEG.csv")

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")
library("dplyr")
library("biomaRt")
library("remotes")
library("nichenetr")

all_TF <- all_DEG
aa <- as.character(all_TF$Gene_symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"GSE12452_all_DEG.csv")


********************************************GSE64634

library("GEOquery")
gse <- getGEO("GSE64634",GSEMatrix = TRUE, destdir = ".",getGPL = T, AnnotGPL = T)
exprs <- exprs(gse[[1]])
fdata <- fData(gse[[1]])
explan <- data.frame(exprs)
explan$ID <- fdata$ID
explan$symbol <- fdata$"Gene symbol"
explan$gene_ID <- fdata$"Gene ID" 
head(explan)
express_2 <- explan
express_2$symbol <- data.frame(sapply(express_2$symbol,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$gene_ID <- data.frame(sapply(express_2$gene_ID,function(x) unlist(strsplit(x,"///"))[1]), stringsAsFactors = F)[,1]
express_2$rowMeans <- apply(express_2[,1:16],1,function(x) mean(as.numeric(x),na.rm = T))
head(express_2)
express_2 <- express_2[order(express_2$rowMeans,decreasing = T),]
express_na <- na.omit(express_2)

explan_final <- express_na[express_na$symbol !="" ,]

explan_final <- express_na[!duplicated(express_na$symbol),]
nrow(explan_final)
colnames(explan_final) <- c("Normal1","Normal2","Normal3","Normal4","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","array_ID","Gene_symbol","gene_ID","rowMeans")
rownames(explan_final) <- explan_final[,18]
count_anno <- explan_final[,1:16]
head(count_anno)
write.csv(count_anno, "GSE64634_count_anno.csv")
count_anno_without_arrayID <- explan_final[,1:19]
head(count_anno_without_arrayID)
write.csv(count_anno_without_arrayID, "GSE64634_count_anno_without_arrayID.csv")

group_list <- factor(c(rep('Normal',4),rep('NPC_tumor',12)), levels=c("Normal","NPC_tumor")) 
table(group_list)

library(limma)
library(data.table)
rt <- count_anno
rt <- as.matrix(rt)

dimnames <- list(rownames(rt),colnames(rt))
data <- matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data <- avereps(data)
data <- as.data.frame(data)
nrow(data)
data <- data[rowMeans(data)>0.5,]

design <- model.matrix(~factor(group_list))
design
fit <- lmFit(data,design)
fit <- eBayes(fit)
DEG <- topTable(fit,coef=2,n=Inf)
head(DEG)

all_DEG <- merge(count_anno_without_arrayID,DEG, by = "row.names", all=FALSE)
head(all_DEG)
write.csv(all_DEG, "GSE64634_all_DEG.csv")

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")
library("dplyr")
library("biomaRt")
library("remotes")
library("nichenetr")

all_TF <- all_DEG
aa <- as.character(all_TF$Gene_symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"GSE64634_all_DEG.csv")


***************************************************************GSE12452 DEG*************************
gsea_data <- read.csv("GSE12452_count_anno.csv")
hg19 <- read.csv("hg19_Ensembl_id.csv")
gsea_data$Description <- c("NULL")

nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
gsea_data <- gsea_data[,c(1,43,2:42)]
colnames(gsea_data) <- c("NAME","Description","Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Normal7","Normal8","Normal9","Normal10","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25","NPC26","NPC27","NPC28","NPC29","NPC30","NPC31")
set <- gsea_data[!duplicated(gsea_data$NAME),]
hg19_data <- merge(gsea_data,hg19,by.x="NAME", by.y = "symbol", all=FALSE)
hg19_data <- hg19_data[,c(44,2:43)]
colnames(hg19_data) <- c("NAME","Description","Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Normal7","Normal8","Normal9","Normal10","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25","NPC26","NPC27","NPC28","NPC29","NPC30","NPC31")

write.csv(hg19_data,file="gsea_data.csv")

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_h -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c1 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c2 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c3 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c4 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c5 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c6 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE12452_tumor_vs_normal_c7 -gui false


***************************************************************GSE13597 DEG*************************
gsea_data <- read.csv("GSE13597_count_anno.csv")
hg19 <- read.csv("./GSE12452/hg19_Ensembl_id.csv")
gsea_data$Description <- c("NULL")

nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
gsea_data <- gsea_data[,c(1,30,2:29)]
colnames(gsea_data) <- c("NAME","Description","Normal1","Normal2","Normal3","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25")
hg19_data <- merge(gsea_data,hg19, by.x="NAME", by.y = "symbol", all=FALSE)
hg19_data <- hg19_data[,c(31,2,6:30,3:5)]
colnames(hg19_data) <- c("NAME","Description","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","NPC13","NPC14","NPC15","NPC16","NPC17","NPC18","NPC19","NPC20","NPC21","NPC22","NPC23","NPC24","NPC25","Normal1","Normal2","Normal3")

write.csv(hg19_data,file="gsea_data.csv")


java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE13597_tumor_vs_normal_h -gui false



***************************************************************GSE64634 DEG*************************
gsea_data <- read.csv("GSE64634_count_anno.csv")
hg19 <- read.csv("./GSE12452/hg19_Ensembl_id.csv")
gsea_data$Description <- c("NULL")

nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
gsea_data <- gsea_data[,c(1,18,2:17)]
colnames(gsea_data) <- c("NAME","Description","Normal1","Normal2","Normal3","Normal4","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12")
hg19_data <- merge(gsea_data,hg19, by.x="NAME", by.y = "symbol", all=FALSE)
hg19_data <- hg19_data[,c(19,2,7:18,3:6)]
head(hg19_data)
colnames(hg19_data) <- c("NAME","Description","NPC1","NPC2","NPC3","NPC4","NPC5","NPC6","NPC7","NPC8","NPC9","NPC10","NPC11","NPC12","Normal1","Normal2","Normal3","Normal4")

write.csv(hg19_data,file="gsea_data.csv")

java -cp ./gsea/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSE64634_tumor_vs_normal_h -gui false



********************************************************primary_TMPL vs primary_TMP******************************
********************************************************primary_TMPL vs primary_TMP******************************
********************************************************primary_TMPL vs primary_TMP******************************
********************************************************primary_TMPL vs primary_TMP******************************

#R3.6

project <- c("primary_LMP1")
deal_design <- c("LMP1","ctrl")
significant_cutoff <- c(0.5)
organism <- "mouse"
file_path <- "."
sample_sampletable.path <- "./star_mm10_out"
load(file="./ebg_mm10.RData")
load(file="./txdb_mm10.RData")  
deal_case <- deal_design[1]
dela_control <- deal_design[2]
deal_names <- paste(deal_case,dela_control,sep="_v_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"00_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"01_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"01_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"01_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"02_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"03",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"04",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"05",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"06",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"07",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"08",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"09",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(DESeq2)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})

if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
	load(file="./ebg_mm10.RData")
	load(file="./txdb_mm10.RData")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
	load(file="./ebg_GRCh38.RData")
	load(file="./txdb_GRCh38.RData")
}

indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sample_sampletable_primary.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",   
singleEnd=FALSE,  
ignore.strand=TRUE,
fragments=FALSE) 

colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)

rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(vsd$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_vsd)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()

#normalize by TPM
exons<-as.data.frame(ebg)
exon_length<-data.frame(gene_name=c(0),length=c(0))
for(i in 1:nrow(assay(se))){
  exonsfromgene<-exons[which(exons$group==i),]
  length<-sum(exonsfromgene$width)
  gene_name<-exonsfromgene[1,2]
  exon_length<-rbind(exon_length,c(gene_name,length))
}
exon_length<-exon_length[-1,]
whole_length<-as.data.frame(round(colSums(assay(dds))/1e6,1))
countM<-as.data.frame(assay(se))
FPKM<-countM
for(i in 1:nrow(FPKM)){
  gene<-rownames(FPKM[i,])
  length<-as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
  for(j in 1:ncol(FPKM)){
    FPKM[i,j]<-FPKM[i,j]/(whole_length[j,]*length*0.001) 
  }
}
fpkmToTpm <- function(FPKM){
    exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- sampleTable$sample
colnames(count) <- sampleTable$sample
count_and_tpm <- cbind(count,tpm)
write.csv(tpm,file=tpm_csv)
#DESeq2  DEG
dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]
#####KEGG 富集
upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)
KEGG_pro_enhan_up <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")
write.csv(KEGG_pro_enhan_up, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_UP_genes.csv")
pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()

downres_1 <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

KEGG_pro_enhan_down <- setReadable(KEGGdownres_1, "org.Mm.eg.db", keyType="ENTREZID")
write.csv(KEGG_pro_enhan_down, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_DOWN_genes.csv")
pdf(file = "primary_LMP1_v_ctrl_KEGG_UP_top15.pdf")
dotplot(KEGGupres_1, showCategory=15)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

#####GO
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=5) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = "primary_LMP1_v_ctrl_GO_DOWN_top5.pdf")
ff1
dev.off()

###GSEA
java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_h -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c1 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c2 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c3 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c4 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c5 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c6 -gui false

java -cp ./gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx ./gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip ./gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c7 -gui false

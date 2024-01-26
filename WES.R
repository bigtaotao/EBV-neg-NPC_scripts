$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1A $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1A $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1A $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#R4.2.2
library(ComplexHeatmap)
library(maftools)
library(RColorBrewer)

#Get the patient's SNVs matrix
neg_CNA <- mcreadRDS("./neg_CNA_pttest_filter.rds",mc.cores=20)
pos_CNA <- mcreadRDS("./pos_CNA_pttest_filter.rds",mc.cores=20)

genes = c("TP53","FOS","CDKN2A","CHD4","CREBBP","MDM2","ATM","CHD3","MTA2",
    "SMAD4","TGFBR2","SMAD2","ACVR1B",
    "NOTCH1","NCOR2","DVL2",
    "PIK3CA","PTEN","PIK3R1","INPPL1",
    "ERBB4","RIT1","FGFR3","KIT","FGFR2",
    "TRAF3","NFKBIA","TRAF2","CYLD","NLRC5",
    "MYC","ZFHX3","BCL9"
    )

ebv_neg_maf_with_cnv <- read.maf(maf="./snp_bld_tnscope.maf",
    cnTable = neg_CNA)
ebv_pos_maf_with_cnv <- read.maf(maf="./snp_hknpc_tnscope.maf",
    cnTable = pos_CNA)

oncoplot(maf = ebv_pos_maf_with_cnv,
    draw_titv = FALSE,  #TiTv plot
    showTitle = TRUE,
    logColBar = FALSE,
    bgCol="#F5F5F5",
    genes=genes,
    keepGeneOrder = TRUE,
    writeMatrix = TRUE
    )

oncoplot(maf = ebv_neg_maf_with_cnv,
    draw_titv = FALSE,  #TiTv plot
    showTitle = TRUE,
    logColBar = FALSE,
    bgCol="#F5F5F5",
    genes=genes,
    keepGeneOrder = TRUE,
    writeMatrix = TRUE
    )

#Get the TMB matrix of patients 
ebv_neg_clinical <- read.csv("./bld_cli.csv",row.names=1)
ebv_pos_clinical <- read.csv("./hknpc_cli.csv",row.names=1)

ebv_neg_maf <- read.maf(maf="./snp_bld_tnscope.maf",
   clinicalData = ebv_neg_clinical,cnLevel="all")
ebv_pos_maf <- read.maf(maf="./snp_hknpc_tnscope.maf",
   clinicalData = ebv_pos_clinical,cnLevel="all")

tmb_ebv_pos <- tmb(maf = ebv_pos_maf,captureSize=1,logScale=FALSE)
tmb_ebv_neg <- tmb(maf = ebv_neg_maf,captureSize=1,logScale=FALSE)

tmb_ebv_pos <- as.data.frame(tmb_ebv_pos)
tmb_ebv_neg <- as.data.frame(tmb_ebv_neg)

ebv_neg_clinical <- merge(ebv_neg_clinical,tmb_ebv_neg[,1:2],by="Tumor_Sample_Barcode")
ebv_pos_clinical <- merge(ebv_pos_clinical,tmb_ebv_pos[,1:2],by="Tumor_Sample_Barcode")

colnames(ebv_neg_clinical)[7] <- "TMB"
colnames(ebv_pos_clinical)[7] <- "TMB"

write.csv(ebv_neg_clinical,"./EBV_neg_clinicla_with_TMB.csv")
write.csv(ebv_pos_clinical,"./EBV_pos_clinicla_with_TMB.csv")

#Draw oncoplot 
ebv_neg_clinical <- read.csv("./EBV_neg_clinicla_with_TMB.csv",row.names=1)
ebv_pos_clinical <- read.csv("./EBV_pos_clinicla_with_TMB.csv",row.names=1)

mat_pos = read.table("./onco_matrix_pos.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mat_pos[is.na(mat_pos)] = ""
mat_pos[1:3, 1:3]

mat_neg = read.table("./onco_matrix_neg.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mat_neg[is.na(mat_neg)] = ""
mat_neg[1:3, 1:3]
mat_neg$bld4_output_tnscope <- ""

ol_patinets_neg <- intersect(colnames(mat_neg),ebv_neg_clinical$Tumor_Sample_Barcode) 
mat_neg <- mat_neg[,ol_patinets_neg]
ebv_neg_clinical <- ebv_neg_clinical[which(ebv_neg_clinical$Tumor_Sample_Barcode %in% ol_patinets_neg),]

ebv_pos_clinical$Tumor_Sample_Barcode <- gsub("-","\\.",ebv_pos_clinical$Tumor_Sample_Barcode)
ol_patinets_pos <- intersect(colnames(mat_pos),ebv_pos_clinical$Tumor_Sample_Barcode) 
mat_pos <- mat_pos[,ol_patinets_pos]
ebv_pos_clinical <- ebv_pos_clinical[which(ebv_pos_clinical$Tumor_Sample_Barcode %in% ol_patinets_pos),]

col = RColorBrewer::brewer.pal(n = 12, name = 'Paired')
names(col) = c("Shallow_Del","Deep_Del",'Multi_Hit','Missense_Mutation','Gain','Amp', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#F5F5F5", col = NA))
    },
    Amp = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Amp"], col = NA))
    },
    Gain = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Gain"], col = NA))
    },
    Deep_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Deep_Del"], col = NA))
    },
    Shallow_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Shallow_Del"], col = NA))
    },
    Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Multi_Hit"], col = NA))
    },
    Missense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Missense_Mutation"], col = NA))
    },
    Nonsense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
    },
    Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
    },
    In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["In_Frame_Ins"], col = NA))
    },
    Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["Splice_Site"], col = NA))
    },
    In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(4, "pt"), h*0.33, 
            gp = gpar(fill = col["In_Frame_Del"], col = NA))
    }
)


#add clinical infomation
library(circlize)
recurrence <- brewer.pal(n = 8, name = 'Pastel1')[1:2]
names(recurrence) <- c("1","0")

smoke=c(brewer.pal(n = 8, name = 'Set3')[c(5,4)],"grey")
smoke=setNames(smoke,c("No","Yes","NIL"))

gender=brewer.pal(n = 8, name = 'Pastel1')[3:4]
names(gender) <- c("male","female")

age=brewer.pal(n = 8, name = 'YlGn')[c(2,4,6)]
age=setNames(age,c("<40","40-60",">60"))

stage=c(brewer.pal(n = 8, name = 'OrRd')[2:8],"grey")
stage=setNames(stage,c("I","II","III","IV","IVA","IVB","IVC","Nil"))

clinical_color <- list(recurrent=recurrence,smoke=smoke,gender=gender,age=age,stage=stage)

ha_H_pos <- HeatmapAnnotation(TMB=anno_barplot(ebv_pos_clinical$TMB,
                                         border=FALSE, gp = gpar(fill="grey")),
                          recurrent=ebv_pos_clinical$recurrent,
                          gender=ebv_pos_clinical$gender,
                          smoke=ebv_pos_clinical$smoke,
                          age=ebv_pos_clinical$age,
                          stage=ebv_pos_clinical$stage,
                          col = clinical_color,
                        show_annotation_name = TRUE,

                        annotation_name_gp = gpar(fontsize = 10))

ha_H_neg <- HeatmapAnnotation(TMB=anno_barplot(ebv_neg_clinical$TMB,
                                         border=FALSE, gp = gpar(fill="grey")),
                          recurrent=ebv_neg_clinical$recurrent,
                          gender=ebv_neg_clinical$gender,
                          smoke=ebv_neg_clinical$smoke,
                          age=ebv_neg_clinical$age,
                          stage=ebv_neg_clinical$stage,
                          col = clinical_color,
                        show_annotation_name = TRUE,
                        annotation_name_gp = gpar(fontsize = 10))


genes = c("TP53","FOS","CDKN2A","CHD4","CREBBP","MDM2","ATM","CHD3","MTA2",
    "SMAD4","TGFBR2","SMAD2","ACVR1B",
    "NOTCH1","NCOR2","DVL2",
    "PIK3CA","PTEN","PIK3R1","INPPL1",
    "ERBB4","RIT1","FGFR3","KIT","FGFR2",
    "TRAF3","NFKBIA","TRAF2","CYLD","NLRC5",
    "MYC","ZFHX3","BCL9"
    )

mat_pos=mat_pos[genes,ebv_pos_clinical$Tumor_Sample_Barcode]
p_pos=oncoPrint(mat_pos,
          top_annotation = ha_H_pos,
          alter_fun = alter_fun, col = col,  
          remove_empty_columns = FALSE, remove_empty_rows = FALSE,
          pct_side = "right", row_names_side = "left",
          row_order = genes,
          column_title = "pos"
)

mat_neg=mat_neg[genes,ebv_neg_clinical$Tumor_Sample_Barcode]
p_neg=oncoPrint(mat_neg,
          top_annotation = ha_H_neg,
          alter_fun = alter_fun, col = col,  
          remove_empty_columns = FALSE, remove_empty_rows = FALSE,
          show_row_names = FALSE,
          pct_side = "right", row_names_side = "left",
          row_order = genes,
          column_title = "neg"
)

ht_list=p_pos+p_neg

pdf(height=10,"./all_oncoplot_v9.pdf")
draw(ht_list, row_km = 2, #row_split = direction,
     column_title = "", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, 
     heatmap_legend_side = "bottom")
dev.off()


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pdf("./TP53_lollipopPlot_pos.pdf")
lollipopPlot(maf = ebv_pos_maf, gene = 'TP53', AACol = 'AAChange.ensGene', showMutationRate = TRUE,labelPos = "all")
dev.off()
pdf("./TP53_lollipopPlot_neg.pdf")
lollipopPlot(maf = ebv_neg_maf, gene = 'TP53', AACol = 'AAChange.ensGene', showMutationRate = TRUE,labelPos = "all")
dev.off()

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig1.C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

ebv_pos <- maftools::readGistic(
    gisticAllLesionsFile = "cnv_hknpc.segments_result/all_lesions.conf_90.txt",
    gisticAmpGenesFile = "cnv_hknpc.segments_result/amp_genes.conf_90.txt",
    gisticDelGenesFile = "cnv_hknpc.segments_result/del_genes.conf_90.txt",
    gisticScoresFile = "cnv_hknpc.segments_result/scores.gistic",
    isTCGA = TRUE,
    verbose = FALSE
)

ebv_neg <- maftools::readGistic(
    gisticAllLesionsFile = "cnv_bld.segments_result/all_lesions.conf_90.txt",
    gisticAmpGenesFile = "cnv_bld.segments_result/amp_genes.conf_90.txt",
    gisticDelGenesFile = "cnv_bld.segments_result/del_genes.conf_90.txt",
    gisticScoresFile = "cnv_bld.segments_result/scores.gistic",
    isTCGA = TRUE,
    verbose = FALSE
)

pdf(width=10,"./gisticChromPlot_pos.pdf")
gisticChromPlot(gistic = ebv_pos, fdrCutOff = 0.1,markBands = "all",cytobandTxtSize=1,mutGenesTxtSize=3,txtSize=0.7)
dev.off()
pdf(width=10,"./gisticChromPlot_neg.pdf")
gisticChromPlot(gistic = ebv_neg, fdrCutOff = 0.1,markBands = "all",cytobandTxtSize=1,mutGenesTxtSize=1,txtSize=0.7)
dev.off()


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS1 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS1 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS1 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

hnsc_maf <- read.maf(maf = "./burden/tcga-hnsc.maf", 
    clinicalData = NULL)

stad_maf <- read.maf(maf = "./burden/tcga-stad.maf", 
    clinicalData = NULL)

tmb_hnsc <- tmb(maf = hnsc_maf,captureSize=1,logScale=FALSE)
tmb_stad <- tmb(maf = stad_maf,captureSize=1,logScale=FALSE)
tmb_ebv_pos <- tmb(maf = ebv_pos_maf,captureSize=1,logScale=FALSE)
tmb_ebv_neg <- tmb(maf = ebv_neg_maf,captureSize=1,logScale=FALSE)

tmb_hnsc <- as.data.frame(tmb_hnsc)
tmb_stad <- as.data.frame(tmb_stad)
tmb_ebv_pos <- as.data.frame(tmb_ebv_pos)
tmb_ebv_neg <- as.data.frame(tmb_ebv_neg)

tmb_hnsc$Tumor_Sample_Barcode <- substr(tmb_hnsc$Tumor_Sample_Barcode,1,15)
tmb_stad$Tumor_Sample_Barcode <- substr(tmb_stad$Tumor_Sample_Barcode,1,15)

stad_ebv_pos_id <- c("TCGA-B7-5818-01","TCGA-BR-4253-01","TCGA-BR-6455-01","TCGA-BR-6706-01","TCGA-BR-6707-01","TCGA-BR-7196-01","TCGA-BR-7958-01","TCGA-BR-8285-01","TCGA-BR-8366-01","TCGA-BR-8381-01","TCGA-BR-8589-01","TCGA-BR-8676-01","TCGA-BR-8686-01","TCGA-BR-A4J4-01","TCGA-CD-5801-01","TCGA-CG-5722-01","TCGA-D7-5577-01","TCGA-D7-8570-01","TCGA-D7-8573-01","TCGA-D7-A4YX-01","TCGA-FP-7916-01","TCGA-FP-7998-01","TCGA-HU-8608-01","TCGA-HU-A4G2-01","TCGA-HU-A4G6-01","TCGA-HU-A4H0-01")

tmb_stad$group <- ifelse(tmb_stad$Tumor_Sample_Barcode %in% stad_ebv_pos_id,"EBV_pos","EBV_neg")
ggboxplot(tmb_stad,x="group",y="total_perMB",fill="group") + 
    stat_compare_means(comparisons=list(c("EBV_pos","EBV_neg")),method="t.test")

hnsc_cli <- read.table("./burden/hnsc-cli.matrix",header=T,sep="\t")
hnsc_ebv_pos_id <- hnsc_cli[which(hnsc_cli$hpv_status_by_ish_testing == "Positive" & hnsc_cli$hpv_status_by_p16_testing == "Positive"),]$sampleID
hnsc_ebv_neg_id <- hnsc_cli[which(hnsc_cli$hpv_status_by_ish_testing == "Negative" & hnsc_cli$hpv_status_by_p16_testing == "Negative"),]$sampleID

tmb_hnsc <- tmb_hnsc[which(tmb_hnsc$Tumor_Sample_Barcode %in% c(hnsc_ebv_pos_id,hnsc_ebv_neg_id)),]
tmb_hnsc$group <- ifelse(tmb_hnsc$Tumor_Sample_Barcode %in% hnsc_ebv_pos_id,"EBV_pos","EBV_neg")
ggboxplot(tmb_hnsc,x="group",y="total_perMB_log",fill="group") + 
    stat_compare_means(comparisons=list(c("EBV_pos","EBV_neg")),method="t.test")

tmb_ebv_pos$group <- "EBV_pos"
tmb_ebv_neg$group <- "EBV_neg"
tmb_ebv <- rbind(tmb_ebv_pos,tmb_ebv_neg)
tmb_ebv$group <- factor(tmb_ebv$group,levels=c("EBV_neg","EBV_pos"))
ggboxplot(tmb_ebv,x="group",y="total_perMB_log",fill="group") + 
    stat_compare_means(comparisons=list(c("EBV_pos","EBV_neg")),method="t.test")

tmb_stad$group <- paste0("STAD_",tmb_stad$group)
tmb_hnsc$group <- paste0("HNSC_",tmb_hnsc$group)
tmb_ebv$group <- paste0("NPC_",tmb_ebv$group)
tmb_all <- rbind(tmb_stad,tmb_hnsc,tmb_ebv)
tmb_all$group <- factor(tmb_all$group,levels=c("NPC_EBV_neg","NPC_EBV_pos","STAD_EBV_neg","STAD_EBV_pos","HNSC_EBV_neg","HNSC_EBV_pos"))

pdf("./TMB_all_cancer.pdf")
ggviolin(tmb_all,x="group",y="total_perMB_log",fill="group",add = "boxplot", add.params = list(fill = "white",alpha=0.5),palette="Paired") + 
    stat_compare_means(comparisons=list(c("NPC_EBV_neg","NPC_EBV_pos"),c("STAD_EBV_neg","STAD_EBV_pos"),c("HNSC_EBV_neg","HNSC_EBV_pos")),method="t.test",label="p.signif")
dev.off()


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2.F $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2.F $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Fig2.F $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


mat_pos=mat_pos[c("TP53","CDKN2A","MYC"),ebv_pos_clinical$Tumor_Sample_Barcode]
p_pos=oncoPrint(mat_pos,
          top_annotation = NULL,
          right_annotation = NULL,
          alter_fun = alter_fun, col = col,  
          remove_empty_columns = FALSE, remove_empty_rows = FALSE,
          pct_side = "right", row_names_side = "left",
          row_order = c("TP53","CDKN2A","MYC"),
          column_title = "pos"
)


mat_neg=mat_neg[c("TP53","CDKN2A","MYC"),ebv_neg_clinical$Tumor_Sample_Barcode]
p_neg=oncoPrint(mat_neg,
          top_annotation = NULL,
          right_annotation = NULL,
          alter_fun = alter_fun, col = col,  
          remove_empty_columns = FALSE, remove_empty_rows = FALSE,
          show_row_names = TRUE,
          pct_side = "right", row_names_side = "left",
          row_order = c("TP53","CDKN2A","MYC"),
          column_title = "neg"
)

pdf(height=4,"./TCM_oncoplot_pos_v4.pdf")
p_pos
dev.off()

pdf(height=4,"./TCM_oncoplot_neg_v4.pdf")
p_neg
dev.off()

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
### add pos neg gc data
if(T){
  rm(list = ls())
  library(data.table)
  library(dplyr)
  library(ggplot2)
  fc_npc = fread("./data/GSE53819.top.table.tsv")
  npc_tr2 = fc_npc %>% filter(Gene.symbol == "TGFBR2")
  npc_logfc = mean(npc_tr2$logFC)
  npc_logfc = data.frame(Tumor_type = "EBV+ NPC",foldchange = npc_logfc)

  # plotdata = rbind(npc_logfc,hnsc_logfc,stad_logfc,gc_pos_logfc,gc_neg_logfc)
  plotdata = npc_logfc
  plotdata$foldchange = 2^plotdata$foldchange
  plotdata$group = "Tumor"
  pn = data.frame(Tumor_type = plotdata$Tumor_type,group = "Normal",foldchange = 1)
  plotdata = rbind(plotdata,pn)
  plotdata$group = factor(plotdata$group,levels = c("Tumor","Normal"))
  p = ggplot(data = plotdata,aes(x = Tumor_type,y = foldchange,fill = group))+
    geom_bar(stat = 'identity',position = "dodge",width=0.4)+
    theme_classic()+
    xlab("Tumor type")+
    ylab("Foldchange")
  ## split pos neg
  p = ggplot(data = plotdata,aes(x = group,y = foldchange,fill = group))+
    geom_bar(stat = 'identity',position = "dodge",width=0.4)+
    theme_classic()+
    xlab("Sample type")+
    ylab("Foldchange")
  ggsave(filename = "./output/21.TGFBR2_foldchange_addGC.pdf",width = 5,height = 5)
}

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.D $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.D $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FigS6.D $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


if(T){
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(paletteer)
  library(data.table)
  library(survminer)
  library(survival)

  get_cond = function(exp_data){
    tgfbr2_exp = c(t(exp_data[row.names(exp_data) == "TGFBR2",]))
    tgfbr2_quantile = quantile(tgfbr2_exp) 
    tgfbr2_low = colnames(exp_data)[tgfbr2_exp < tgfbr2_quantile[2]]
    tgfbr2_high = colnames(exp_data)[tgfbr2_exp > tgfbr2_quantile[4]]
    tgfbr2_low = data.frame(sample = tgfbr2_low,TGFBR2 = "low")
    tgfbr2_high = data.frame(sample = tgfbr2_high,TGFBR2= "high")
    tgfbr2_mid = colnames(exp_data)[tgfbr2_exp <= tgfbr2_quantile[4] & tgfbr2_exp >= tgfbr2_quantile[2]]
    tgfbr2_mid = data.frame(sample = tgfbr2_mid,TGFBR2 = "mid")
    cond = rbind(tgfbr2_low,tgfbr2_high,tgfbr2_mid)
    # rownames(cond) = cond$sample
    # cond$group = cond$condition
    # cond = cond[,-1]
    return(cond)
  }
  #EBV+ NPC
  exp_data = read.csv("./data/ebv_NPC_genes.fpkm_table.csv",row.names = 1)
  os_data = read.csv("./data/PFS.csv") %>% dplyr::select("case","status","time")
  colnames(exp_data) = str_sub(colnames(exp_data),end = -3)
  exp_data = exp_data[,colnames(exp_data) %in% os_data$case]
  
  condition = get_cond(exp_data)
  survdata = left_join(condition,os_data,by = c("sample" = "case"))
  fit_km =surv_fit(Surv(time,status) ~TGFBR2, data = survdata)
  pdf(paste0("./output/19.NPC_TGFBR2_survival.pdf"),width = 5,height = 5)
  print(ggsurvplot(fit_km, data = survdata,surv.median.line = "hv",  # 增加中位生存时间
                   conf.int = F,risk.table = F,title = "TGFBR2 high-low survival status", pval = TRUE))
  dev.off()
}

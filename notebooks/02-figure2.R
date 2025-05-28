# 1. LOAD LIBRARIES
library(Seurat)          # FindMarkers, FetchData
library(dplyr)           # data wrangling
library(ggplot2)         # custom plotting
library(ComplexHeatmap)  # module heatmaps
library(circlize)        # color mapping for heatmaps
library(UCell)           # AddModuleScore_UCell
library(ggpubr)          # stat_compare_means

# 2. PATHS & I/O
data_file    <- "data/HCCprocesseddata.rds"
deg_dir      <- "data/DEGs"                 
results_dir  <- "results"
dir.create(results_dir, showWarnings = FALSE)

# 3. READ IN OBJECT + DEGs
se  <- readRDS(data_file)
# list of cell types
celltypes <- c("Monocyte","CD8+ T","CD4+ T","NK","B")

# load each DEG table into a single data.frame
deg_list <- lapply(celltypes, function(ct) {
  read.csv(file.path(deg_dir, paste0("MASTPositive_DEGs_", ct, ".csv"))) %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
    mutate(
      Celltype  = ct,
      Direction = if_else(avg_log2FC > 0, "Up-regulated","Down-regulated")
    )
})
degs <- bind_rows(deg_list)

# 4. FIG 2A: BAR CHART OF #DEGs
degs_summary <- degs %>%
  count(Celltype, Direction) %>%
  group_by(Celltype) %>%
  mutate(Total = sum(n))

p2A <- ggplot(degs_summary, aes(x = Celltype, y = n, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c("Down-regulated"="#CFAFD4","Up-regulated"="#74AED4")) +
  labs(y="Nr of expression differences") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggsave(file.path(results_dir,"Fig2A_DEG_counts.png"), p2A, width=6, height=4)

# 5. FIG 2B: CLUSTERED HEATMAP OF MODULES
# build a wide matrix of avg_log2FC for each gene×celltype
heatmat <- degs %>%
  select(Gene=gene, Celltype, avg_log2FC) %>%
  pivot_wider(names_from=Celltype, values_from=avg_log2FC, values_fill=0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# run k-means to get 6 modules
set.seed(42)
km <- kmeans(heatmat, centers=6, nstart=20)
heatmat_mod <- split(as.data.frame(heatmat), km$cluster)

# define a blue↔white↔red color map
col_fun <- colorRamp2(c(-1,0,1), c("#CFAFD4","white","#74AED4"))

# draw each cluster as a separate heatmap block
ht_list <- lapply(names(heatmat_mod), function(mod){
  Heatmap(
    as.matrix(heatmat_mod[[mod]]),
    name = paste0("Module", mod),
    col  = col_fun,
    show_row_names    = FALSE,
    show_column_names = TRUE,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE
  )
})
pdf(file.path(results_dir,"Fig2B_Modules_heatmap.pdf"), width=5, height=8)
draw(Reduce(`%v%`, ht_list), heatmap_legend_side="bottom")
dev.off()

# 6. PREPARE FOR MODULE VIOLIN PLOTS (C–F)
# define your gene sets (example)
modules <- list(
  UBIQUITIN  = c("UBA2,SAE1,HUWE1,STUB1,UBE4B,ANAPC10,PIAS3,UBE2E3,WWP1,WWP2,UBE2C,ERCC8,UBE2J2,UBE2QL1,UBE2F,UBE2U,DDB1,DDB2,UBOX5,TRIM32,RHOBTB2,FBXW11,MGRN1,NEDD4L,CBLC,PPIL2,CDC26,ANAPC13,RCHY1,HERC4,FBXO2,FBXW8,FBXO4,UBE2S,PRPF19,ANAPC2,ANAPC4,UBE2K,AIRE,BIRC2,BIRC3,XIAP,NHLRC1,UBE2NL,MDM2,MAP3K1,MID1,TRIM37,NEDD4,PRKN,FZR1,UBR5,ANAPC5,ANAPC7,UBE2J1,ANAPC11,PIAS4,UBE2D4,PML,UBE2R2,DET1,FANCL,UBA6,UBE2W,FBXW7,UBE2Q1,KLHL9,SMURF1,BIRC6,UBE2O,COP1,ANAPC1,SMURF2,SIAH1,SKP1,SKP2,,,UBE2Z,,BRCA1,ELOC,ELOB,TRAF6,SKP1P2,UBA1,UBA7,UBE2A,UBE2B,UBE2D1,UBE2D2,UBE2D3,UBE2E1,UBE2E2,UBE2G1,UBE2G2,UBE2H,UBE2I,UBE2L3,UBE2N,UBE3A,VHL,CUL5,ITCH,SYVN1,CUL4B,CUL4A,CUL3,CUL2,CUL1,PIAS1,SOCS1,CBL,CBLB,CDC23,CDC16,HERC3,HERC2,HERC1,BTRC,UBE3B,SOCS3,KLHL13,UBA3,UBE2M,PIAS2,UBE2L6,UBE2Q2,TRIP12,UBE4A,RNF7,UBE3C,KEAP1,CUL7,CDC20,CDC27,CDC34,RBX1"),  # fill in your Module 2 genes
  TGF_BETA   = c("ACVR1,APC,ARID4B,BCAR3,BMP2,BMPR1A,BMPR2,CDH1,CDK9,CDKN1C,CTNNB1,ENG,FKBP1A,FNTA,FURIN,HDAC1,HIPK2,ID1,ID2,ID3,IFNGR2,JUNB,KLF10,LEFTY2,LTBP2,MAP3K7,NCOR2,NOG,PMEPA1,PPM1A,PPP1CA,PPP1R15A,RAB31,RHOA,SERPINE1,SKI,SKIL,SLC20A1,SMAD1,SMAD3,SMAD6,SMAD7,SMURF1,SMURF2,SPTBN1,TGFB1,TGFBR1,TGIF1,THBS1,TJP1,TRIM33,UBE2D3,WWTR1,XIAP"),  # Module 6
  IMMUNITY1  = c("C17orf99,PARP3,IGLL5,KLRF2,MICA,KLRC4-KLRK1,RASGRP1,CD96,IGHV2-70D,IGHV1-69D,IGHV3-64D,,TCIRG1,TUBB4B,BTN3A3,SPON2,MAD2L2,BATF,CEBPG,CD226,MASP2,CTSC,TRAF3IP2,CPLX2,LILRB1,ARID5A,HCST,MALT1,LILRB4,RIPK3,BTN3A2,CD160,CHGA,CORO1A,SHLD3,TREX1,SCN11A,LYST,CD300A,TUSC2,NLRP3,SLAMF6,RASGRP4,CLNK,SH2D1B,MRGPRX2,CLC,EXOSC6,CLU,SLC15A4,CCR6,RNF19B,RAET1E,ADORA2B,CR1,CR1L,CR2,CRK,CSF2RB,PIK3R6,TICAM1,IL23R,SHLD1,CTSG,CTSH,CX3CR1,RAET1L,CD55,DAO,DBH,DENND1B,ACE,AP1G1,DDX1,RNF168,DHX36,NLRP6,AGER,DNASE1,DNASE1L3,JAG1,AHR,ELANE,APLF,UNC13D,EMP2,GAPT,TUBB,ERCC1,F2,F2RL1,FCER1A,FCER1G,FCER2,FCGR1A,FCGR1BP,FCGR2A,FCGR2B,FCGR3A,FCGR3B,FES,FGR,ZBTB1,KLRK1,FOXF1,PAXIP1,FOXJ1,SWAP70,RFTN1,PLEKHM2,UFL1,CLCF1,RIGI,CADM1,FUT7,NCR3,IL4I1,STAP1,GATA1,GATA2,GATA3,GFER,LAT,RABGEF1,GNL1,IGHV8-51-1,IGHV7-81,IGHV6-1,IGHV5-10-1,IGHV5-51,,IGHV4-61,IGHV4-59,IGHV4-39,IGHV4-34,IGHV4-31,,IGHV4-28,IGHV4-4,MILR1,,IGHV3-74,IGHV3-73,IGHV3-72,IGHV3-66,IGHV3-64,IGHV3-53,IGHV3-49,IGHV3-48,IGHV3-43,IGHV3-38,IGHV3-35,IGHV3-33,IGHV3-30,IGHV3-23,IGHV3-21,IGHV3-20,IGHV3-16,IGHV3-15,IGHV3-13,IGHV3-11,,IGHV3-7,IGHV2-70,IGHV2-26,IGHV2-5,IGHV1-69-2,,IGHV1-69,IGHV1-58,IGHV1-45,IGHV1-24,IGHV1-18,,IGHV1-3,IGLC7,GRB2,CD274,GRP,MSH6,TBX21,GZMB,GZMM,ANXA3,NCKAP1L,HFE,ADGRE2,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-E,HLA-F,HLA-G,HLA-H,MR1,HMGB1,HPRT1,AIRE,HPX,HSPA8,HSPD1,ICAM1,CLEC4G,CFI,IFNA2,IFNB1,IGHA1,IGHA2,IGHD,IGHE,IGHG1,IGHG2,IGHG3,IGHG4,IGHM,IGKC,RAET1G,IGLC1,IGLC2,IGLC3,IGLC6,IGLL1,IL1B,IL1R1,IL2,IL2RB,IL4,IL4R,IL6,IL7R,IL9,IL9R,IL10,IL12A,IL12B,IL12RB1,IL13,IL13RA2,IL18,INPP5D,IRF7,ITGAM,ITGB2,JAK3,NCR3LG1,KIF5B,KIR2DL4,KIR3DL1,KIT,KLRB1,KLRC1,KLRC2,KLRC3,KLRD1,ARG1,GPR15LG,CLEC2A,CLEC12B,SCIMP,LAG3,LAMP1,LEP,LGALS9,LIG4,RAB44,LTA,LYN,SH2D1A,ARRB2,SMAD7,MBL2,CD46,MICB,MLH1,MPL,MSH2,MYD88,NBN,NKG7,NOS2,P2RX7,PLA2G3,PRDX1,IL21R,FOXP3,EXOSC3,KMT5B,IRAK4,C1RL,TLR8,PDCD1,IL23A,CYRIB,PDPK1,SERPINB9,PIK3CD,PIK3CG,PIK3R1,PLA2G1B,PLCG2,IL20RB,PMS2,TLR9,TREM2,TREM1,SASH3,SHLD2,TRPM4,SUSD4,RIF1,ARL8B,PPP3CB,PRF1,PRKCD,PRKCZ,KMT2E,CRTAM,AZU1,B2M,SPHK2,DUSP22,HMCES,C12orf4,CD177,PTAFR,IGHV7-4-1,PTGDR,,PTGDS,AICDA,MAVS,PTPN6,SLAMF7,PTPRC,PVR,NECTIN2,SNX6,RAB27A,RAC2,IL21,BCL3,BCL6,BCR,SCART1,S100A13,SERPINB4,SCNN1B,CEACAM1,CCL3,CXCL6,XCL1,MYO1G,NOD2,CARD9,NFKBIZ,CLEC7A,SLAMF1,SLC18A2,SPI1,SPN,STAT5B,STAT6,STX4,STXBP1,STXBP2,STXBP3,SUPT6H,BST2,VAMP2,VAMP7,SYK,ADAM17,MAP3K7,TAP2,BTK,,,TFRC,TGFB1,C1QBP,TLR3,TLR4,SERPING1,C1QA,TNF,C1QB,TNFRSF1B,C1QC,C1R,TP53BP1,C1S,C2,C3,TRAF2,TRAF6,C4A,C4B,C4BPA,C4BPB,C5,C6,TNFSF4,CCR2,C7,TYROBP,C8A,C8B,C8G,C9,UNG,VAV1,WAS,LAT2,NSD2,ZP3,FZD5,ULBP3,ATAD5,SVEP1,NR4A3,ULBP2,ULBP1,NDFIP1,CAMK4,FBXO38,UNC93B1,KDM5D,KLRC4,PRAM1,STX7,SLA2,JAGN1,SANBR,KMT5C,HAVCR2,NDST2,CBL,VAMP8,SNX4,TNFSF13,FADD,SNAP23,IL18RAP,IL18R1,CD84,BCL10,,SLAMF9,RNF8,CD1A,CD1B,FCGR2C,CD1C,CD1D,CD1E,CD2,RSAD2,EXO1,EBAG9,DDX21,FCMR,CD8A,CD19,VAMP3,CD27,SLC22A13,CD28,CD80,NCR1,IL27RA,CD40,CD40LG,CD70,CD74,CD81,GAB2,WDR1"),  # Module 1
  ISG        = c('BST2','CASP1','CMPK2','DDX60','EIF2AK2','EPSTI1','TENT5A','HERC6','IFI35','IFI44','IFI44L',
             'IFIT3','IFITM2','IFITM3','IRF2','IRF7','IRF9','ISG15','ISG20','LAP3','LPAR6','LY6E','MX1',
             'OAS1','OASL','PARP9','PLSCR1','PNPT1','HELZ2','PSMB9','PSME1','RSAD2','SAMD9','SAMD9L',
             'SELL','SP110','STAT2','TRIM14','TRIM25','TRIM5','UBE2L6','APOL6','CD38','CIITA',
             'JAK2','MT2A','MX2','OAS2','OAS3','PSMB10','PTPN6','RBCK1','RNF213','STAT1',
             'TNFSF10','XAF1','ZBP1')     # Modules 4+5
)

# add UCell scores
for (mod in names(modules)) {
  se <- AddModuleScore_UCell(se, features=list(modules[[mod]]), name=mod)
}

# fetch data for plotting
plotdat <- FetchData(se, vars=c(names(modules), "CELLtype","time")) %>%
  mutate(time = factor(recode(time,"Control"="hc"), levels=c("pre","post","hc")))

# 7. FIG 2C–F: VIOLIN PLOTS FOR EACH MODULE
comparison_list <- list(c("pre","post"), c("pre","hc"), c("post","hc"))

for (mod in names(modules)) {
  p <- ggplot(plotdat, aes_string(x="time", y=paste0(mod,"_UCell"))) +
    geom_violin(trim=FALSE, scale="width") +
    geom_jitter(size=0.1, alpha=0.3) +
    stat_compare_means(comparisons=comparison_list, label="p.format") +
    facet_wrap(~ CELLtype, ncol=5) +
    labs(
      title = paste0(mod," score"),
      y     = "AUC score",
      x     = NULL
    ) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1))
  
  ggsave(
    filename = file.path(results_dir, paste0("Fig2_",mod,"_vln.png")),
    plot     = p,
    width    = 8, height = 4
  )
}

# 8. FIG 2G–H: *CellChat* and chord diagrams
# (see your 03-cellchat.R for these panels)

message("Figure 2A–F complete. See results/*.")  

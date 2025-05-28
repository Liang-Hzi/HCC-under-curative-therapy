# ────────────────────────────────────────────────────────────
# 04-figure4.R: Tissue comparison & module scores (Fig 4A–E)
# ────────────────────────────────────────────────────────────

# 1. LIBRARIES
library(Seurat)    # for DimPlot & FetchData
library(dplyr)     # data wrangling
library(ggplot2)   # plotting
library(ggrepel)   # label volcano
library(ggpubr)    # stat_compare_means
library(UCell)     # AddModuleScore_UCell

# 2. I/O
tissue_rds   <- "data/HCCtissue_immunecells.rds"
shared_csv   <- "data/DEGS/monocyte_shared.csv"
results_dir  <- "results"
dir.create(results_dir, showWarnings = FALSE)

# 3. LOAD & PREPARE
tiss <- readRDS(tissue_rds)
# drop two patients
tiss <- subset(tiss, patient %nin% c("HCC01","HCC02"))
tiss$site <- factor(tiss$site, levels=c("Adjacent liver","Tumor"))

# set celltype colors
ct_cols <- c(
  "Hepatocyte"   = "#5A7E4C",
  "Endothelial"  = "#B0D8B7",
  "Fibroblast"   = "#AA7A38",
  "CD8+ T"       = "#ECA8A9",
  "CD4+ T"       = "#D3E2B7",
  "B cell"       = "#74AED4",
  "NK"           = "#F7C97E",
  "Monocyte"     = "#CFAFD4",
  "DC"           = "#8680A6"
)

# ────────────────────────────────────────────────────────────
# Fig 4A: UMAP split by site
# ────────────────────────────────────────────────────────────
p4A <- DimPlot(
  tiss, reduction="umap", split.by="site", group.by="CELLtype",
  cols=ct_cols, pt.size=0.5, label=FALSE
) + theme_minimal() +
    theme(legend.text=element_text(size=9), axis.text=element_text(size=10))
ggsave(file.path(results_dir,"Fig4A_UMAP_site.png"),
       p4A, width=10, height=4)

# ────────────────────────────────────────────────────────────
# Fig 4B: volcano of shared monocyte DEGs
# ────────────────────────────────────────────────────────────
shared <- read.csv(shared_csv)
shared$gene <- shared$X
# assign pathways
IFN     <- strsplit("MT2A/NAMPT/FGL2/IFITM2/…/SAMD9","/")[[1]]
Inflam  <- strsplit("CD55/CYBB/…/IFIH1","/")[[1]]
Macro   <- strsplit("GAPDH/SQSTM1/…/TRIM13","/")[[1]]
OxPhos  <- strsplit("ATP6V1F/…/ATP6V1A","/")[[1]]
shared$pathway <- "Other"
shared$pathway[shared$gene %in% IFN]    <- "Interferon response"
shared$pathway[shared$gene %in% Inflam] <- "Positive regulation of\ncytokine production"
shared$pathway[shared$gene %in% Macro]  <- "Macroautophagy"
shared$pathway[shared$gene %in% OxPhos] <- "Oxidative phosphorylation"

# avoid zeros
shared$p_val_adj[shared$p_val_adj==0] <- 
  min(shared$p_val_adj[shared$p_val_adj>0])*1e-1

p4B <- ggplot(shared, aes(x=avg_log2FC,y=-log10(p_val_adj), color=pathway)) +
  geom_point(alpha=0.7, size=1.3) +
  scale_color_manual(values=c(
    "Interferon response"="#D73027",
    "Positive regulation of\ncytokine production"="#762A83",
    "Macroautophagy"="#4575B4",
    "Oxidative phosphorylation"="#1B9E77",
    "Other"="grey60"
  )) +
  geom_vline(xintercept=c(-0.25,0.25), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(
    data = subset(shared, pathway!="Other"),
    aes(label=gene), size=2.5, max.overlaps=12
  ) +
  labs(
    title="Monocyte: Adjacent liver vs Tumor",
    x="Log2 Fold Change", y="-Log10 FDR", color="Pathway"
  ) +
  theme_minimal() +
  theme(
    legend.position="right",
    text=element_text(size=11)
  )
ggsave(file.path(results_dir,"Fig4B_Volcano_shared.png"),
       p4B, width=6, height=5)

# ────────────────────────────────────────────────────────────
# Fig 4C–E: module scores by site
# ────────────────────────────────────────────────────────────
# define feature sets
LEUK <- strsplit("…/…/…", "/")[[1]]   # your Leukocyte list
IFN  <- strsplit("…/…/…", "/")[[1]]   # your IFN list
UBQ  <- strsplit("…/…/…", "/")[[1]]   # your Ubiquitin list

# compute scores
tiss <- AddModuleScore_UCell(tiss, list(LEUK), name="Leukocyte")
tiss <- AddModuleScore_UCell(tiss, list(IFN),    name="IFNresp")
tiss <- AddModuleScore_UCell(tiss, list(UBQ),    name="Ubq")

df <- FetchData(
  tiss, 
  vars=c("Leukocyte_UCell","IFNresp_UCell","Ubq_UCell","CELLtype","site")
)

# helper to plot
make_vln <- function(var,title){
  ggplot(df, aes_string(x="site",y=paste0(var,"_UCell"),fill="site")) +
    geom_violin(trim=FALSE, scale="width") +
    geom_jitter(size=0.3, alpha=0.4, width=0.2) +
    stat_compare_means(method="wilcox.test", label="p.format") +
    facet_wrap(~CELLtype, ncol=5) +
    scale_fill_manual(values=c("Adjacent liver"="#825CA6","Tumor"="#94221F")) +
    labs(title=title,y="AUC score",x=NULL) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none")
}

# Fig 4C
p4C <- make_vln("Leukocyte","Leukocyte mediated immunity")
ggsave(file.path(results_dir,"Fig4C_Leukocyte.png"),p4C, width=8, height=4)

# Fig 4D
p4D <- make_vln("IFNresp","Interferon response")
ggsave(file.path(results_dir,"Fig4D_IFNresp.png"),p4D, width=8, height=4)

# Fig 4E (only CD4+ T)
sub <- df %>% filter(CELLtype=="CD4+ T")
p4E <- ggplot(sub, aes(x=site,y=Ubq_UCell,fill=site)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_jitter(size=0.3, alpha=0.4, width=0.2) +
  stat_compare_means(method="wilcox.test", label="p.format") +
  labs(title="Ubiquitin mediated proteolysis",y="AUC score",x=NULL) +
  scale_fill_manual(values=c("Adjacent liver"="#825CA6","Tumor"="#94221F")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none")
ggsave(file.path(results_dir,"Fig4E_Ubiquitin.png"),p4E, width=4, height=4)

message("✅ Figure 4A–E complete; see results/*.")  

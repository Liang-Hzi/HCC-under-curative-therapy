# ────────────────────────────────────────────────────────────
# 06-figure6.R: Bulk RNA volcano, TCGA survival, ABCA1 grouping (Fig 6A–E)
# ────────────────────────────────────────────────────────────

# 1. LIBRARIES
library(ggplot2)       # plotting
library(ggrepel)       # volcano labels
library(DESeq2)        # bulk differential expression
library(clusterProfiler) # optional enrichment later
library(survival)      # KM curves
library(survminer)     # ggsurvplot
library(dplyr)         # data wrangling
library(openxlsx)      # read.xlsx
library(ggpubr)        # stat_compare_means

# 2. I/O
bulk_xlsx    <- "data/Mono3markers_23.xlsx" 
clin_csv     <- "data/TCGAphenotype.csv"
expr_csv     <- "data/TCGAmatrix.csv"
seurat_rds   <- "data/HCCprocesseddata.rds"       # for ABCA1 grouping
results_dir  <- "results"
dir.create(results_dir, showWarnings = FALSE)

# ────────────────────────────────────────────────────────────
# Fig 6A: volcano of bulk CD14+ mono3 signature (HCC vs hc)
# ────────────────────────────────────────────────────────────
volc_df <- read.xlsx(bulk_xlsx, 1) %>%
  rename(
    gene        = X,
    log2FC      = log2FoldChange,
    padj        = padj
  ) %>%
  mutate(
    sig = ifelse(abs(log2FC) > 0.25 & padj < 0.05, "highlight","other")
  )

p6A <- ggplot(volc_df, aes(x=log2FC,y=-log10(padj),color=sig)) +
  geom_point(size=2,alpha=0.8) +
  scale_color_manual(values=c("other"="grey80","highlight"="#AA7A38")) +
  geom_vline(xintercept=c(-0.25,0.25), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(data=filter(volc_df,sig=="highlight"),aes(label=gene),size=3) +
  labs(
    title="CD14+ mono3 signature (Bulk RNA-seq)",
    x="log₂ Fold Change (HCC vs hc)",
    y="-log₁₀ FDR"
  ) +
  theme_minimal() + theme(legend.position="none")

ggsave(file.path(results_dir,"Fig6A_bulk_volcano.png"),p6A,width=5,height=4)

# ────────────────────────────────────────────────────────────
# Fig 6B: Kaplan–Meier for TCGA LIHC (median split)
# ────────────────────────────────────────────────────────────
clin <- read.csv(clin_csv, stringsAsFactors=FALSE)
mat  <- read.csv(expr_csv, row.names=1, check.names=FALSE)

# define genes used in signature:
sig_genes <- volc_df$gene[volc_df$sig=="highlight"]
# compute mean signature per sample, correct by CD45:
cd45    <- "PTPRC"
sig_raw <- colMeans(mat[sig_genes, ], na.rm=TRUE)
sig_adj <- sig_raw / mat[cd45, ]

clin$signature <- sig_adj[match(clin$barcode, names(sig_adj))]
clin$group     <- ifelse(clin$signature > median(clin$signature, na.rm=TRUE),
                         "High(n=185)","Low(n=185)")

# fit & plot
fit6B <- survfit(Surv(time/30, status) ~ group, data=clin)
p6B <- ggsurvplot(
  fit6B, data=clin,
  pval=TRUE, risk.table=FALSE,
  palette=c("High(n=185)"="red","Low(n=185)"="blue"),
  legend.title="Group",
  ggtheme=theme_classic(),
  title="Overall Survival (TCGA LIHC)\nCD14+ mono3 signature"
)

ggsave(file.path(results_dir,"Fig6B_TCGA_KM.png"),
       p6B$plot, width=5, height=4)

# ────────────────────────────────────────────────────────────
# Fig 6C–E: ABCA1 high vs low in single-cell
# ────────────────────────────────────────────────────────────
# load and subset to monocytes
se <- readRDS(seurat_rds)
mono_ct <- c("CD14+ mono1","CD14+ mono2","CD14+ mono3","CD14+ mono4","CD16+ mono")
Mono <- subset(se, subset=CELLtype %in% mono_ct)
Mono$time <- factor(recode(Mono$time,"Control"="hc"),
                    levels=c("pre","post","hc"))

# extract ABCA1 expression, define groups
Mono$ABCA1 <- FetchData(Mono, vars="ABCA1")[,1]
Mono$ABCA1_grp <- ifelse(Mono$ABCA1>0,"ABCA1+","ABCA1-")

# prepare df for plotting
dfp <- Mono@meta.data %>%
  count(time, ABCA1_grp) %>%
  group_by(time) %>%
  mutate(pct=n/sum(n)*100) %>%
  ungroup()

# C: box + lines
p6C <- ggplot(dfp, aes(x=time,y=pct,fill=time)) +
  geom_boxplot(alpha=0.6,outlier.shape=NA) +
  geom_jitter(aes(color=ABCA1_grp),size=1,position=position_jitter(width=0.2)) +
  stat_compare_means(aes(label=..p.signif..), method="wilcox.test",
                     comparisons=list(c("pre","post")), label.y=125) +
  scale_fill_manual(values=c("pre"="#F2B342","post"="#5AAA46","hc"="#4F63B0")) +
  labs(title="ABCA1+ vs ABCA1− cell proportions", y="Percentage (%)", x="Time") +
  theme_bw()

ggsave(file.path(results_dir,"Fig6C_ABCA1_props.png"),p6C,width=5,height=4)

# D: bar chart of #DEGs for ABCA1+ vs ABCA1− at each time
# you can compute these directly via FindMarkers, but here’s a manual table:
deg_counts <- data.frame(
  time=c("pre","post","hc"),
  ABCA1.pos=c(612,27,2),
  ABCA1.neg=c(115,9,1)
) %>%
  pivot_longer(-time,names_to="group",values_to="Nr_DEGs")

p6D <- ggplot(deg_counts, aes(x=time,y=Nr_DEGs,fill=group)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("ABCA1.pos"="red","ABCA1.neg"="blue"),
                    labels=c("ABCA1+","ABCA1−")) +
  labs(title="Nr of expression differences", x=NULL, y="Nr of DEGs", fill="Group") +
  theme_classic()

ggsave(file.path(results_dir,"Fig6D_DEG_counts.png"),p6D,width=5,height=4)

# E: volcano for pre-therapy ABCA1− vs ABCA1+
mamark <- FindMarkers(
  subset(Mono, subset=time=="pre"), 
  ident.1="ABCA1-", ident.2="ABCA1+", 
  test.use="MAST", logfc.threshold=0, min.pct=0.1
)
mamark$gene <- rownames(mamark)
mamark$group <- with(mamark,
  ifelse(avg_log2FC>0.25 & p_val_adj<0.05,"ABCA1− up",
  ifelse(avg_log2FC< -0.25 & p_val_adj<0.05,"ABCA1+ up","other"))
)

# pick genes to label
labels <- c("ABCA1","RNF24","PIK3R5","RAB8B","MEF2A","RICTOR")

p6E <- ggplot(mamark, aes(x=avg_log2FC,y=-log10(p_val_adj),color=group)) +
  geom_point(alpha=0.7,size=1.5) +
  scale_color_manual(values=c("other"="grey80", "ABCA1− up"="#C51B8A","ABCA1+ up"="#74AED4")) +
  geom_vline(xintercept=c(-0.25,0.25), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(data=subset(mamark,gene%in%labels),aes(label=gene),size=3) +
  labs(title="ABCA1− vs ABCA1+ (pre-therapy)", x="Log₂ FC", y="-Log₁₀ FDR") +
  theme_minimal() + theme(legend.position="none")

ggsave(file.path(results_dir,"Fig6E_ABCA1_volcano.png"),p6E,width=5,height=4)

message("✅ Figure 6 panels A–E saved under results/")  

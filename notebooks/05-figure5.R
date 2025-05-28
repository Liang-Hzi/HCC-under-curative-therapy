# ────────────────────────────────────────────────────────────
# 05-figure5.R: Milo DA + tissue‐mac subtypes (Fig 5A–G)
# ────────────────────────────────────────────────────────────

# 1. LIBRARIES
library(Seurat)      # DimPlot, FetchData
library(miloR)       # Milo DA
library(ggplot2)     # plotting
library(ggrepel)     # label points
library(dplyr)       # data wrangling
library(UCell)       # module scores
library(ggpubr)      # stat_compare_means

# 2. I/O
mono_rds    <- "data/HCCtissue_immunecells_monocyte.rds"
milo_rds    <- "data/monocyte_milo.rds"            # your saved Milo object
results_dir <- "results"
dir.create(results_dir, showWarnings=FALSE)

# ────────────────────────────────────────────────────────────
# Fig 5A: Milo DA pre vs post
# ────────────────────────────────────────────────────────────
# load a pre-built Milo object (with cells, graph + counts)
milo_obj <- readRDS(milo_rds)

# If you haven’t already:
# milo_obj <- calcNhoods(milo_obj, prop=0.1, k=20, d=50)
# milo_obj <- buildNhoodGraph(milo_obj)

# prepare a design data.frame from your cell‐metadata
colData(milo_obj)$time <- factor(colData(milo_obj)$time,
                                 levels=c("pre","post","hc"))
design <- data.frame(time=colData(milo_obj)$time,
                     row.names=colnames(milo_obj))

# test DA (pre vs post)
da <- testNhoods(
  milo_obj,
  design=design,
  design_formula=~time,
  contrast=c("time","post","pre"),
  reduced=FALSE
)

# extract results
da_dt <- as.data.frame(da@metadata$stats) %>%
  rename(logFC = logFC, PValue = p.value, FDR = adj.p.value)

# volcano‐style DA plot
p5A <- ggplot(da_dt, aes(x=logFC, y=-log10(PValue), color=FDR<0.1)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("grey70","tomato")) +
  labs(title="DA: post vs pre (Milo)", x="Log₂FC", y="-Log₁₀ P") +
  theme_minimal()
ggsave(file.path(results_dir,"Fig5A_Milo_DA.png"),
       p5A, width=5, height=4)

# ────────────────────────────────────────────────────────────
# Fig 5B: UMAP of monocyte subtypes in tissue
# ────────────────────────────────────────────────────────────
mono <- readRDS(mono_rds)
# subset to only your tissue‐monocyte object
Idents(mono) <- "subtype"
cols <- c(
  "Mac1"="#B39DDB","Mac2"="#9575CD","Mac3"="#BA68C8",
  "Mac4"="#7E57C2","Mac5"="#AB47BC","Mono1"="#EC407A","Mono2"="#CFAFD4"
)
p5B <- DimPlot(
  mono, reduction="umap", group.by="subtype",
  cols=cols, label=FALSE, pt.size=0.5
) + theme_minimal() +
    theme(legend.text=element_text(size=9),
          axis.text=element_text(size=10))
ggsave(file.path(results_dir,"Fig5B_UMAP_tissue_mono.png"),
       p5B, width=5, height=4)

# ────────────────────────────────────────────────────────────
# Fig 5C–E–F–G: violin plots for your signatures
# (CD14+ mono3 signature, M2, IL-10)
# ────────────────────────────────────────────────────────────
# replace these with your real gene lists
mono3_genes <- scan("data/Mono3markers_94.csv", what=character(), sep=",")
M2_genes    <- c("ARG1","CD163","IL10","…")     # your M2 list
IL10_genes  <- c("TNF","CCL2","IL10","…")       # your IL-10 list

# compute scores
mono <- AddModuleScore_UCell(mono, list(mono3_genes), name="C3")
mono <- AddModuleScore_UCell(mono, list(M2_genes),    name="M2")
mono <- AddModuleScore_UCell(mono, list(IL10_genes),  name="IL10")

df <- FetchData(
  mono,
  vars=c("C3_UCell","M2_UCell","IL10_UCell","subtype")
)

plot_violin <- function(df, var, title){
  ggplot(df, aes_string(x="subtype", y=paste0(var,"_UCell"), fill="subtype")) +
    geom_violin(trim=FALSE, scale="width") +
    geom_jitter(size=0.1, alpha=0.3, width=0.2) +
    stat_compare_means(comparisons=list(
      c("Mac1","Mac5"),c("Mac2","Mac5"),c("Mac3","Mac5"),
      c("Mac4","Mac5"),c("Mono1","Mac5"),c("Mono2","Mac5")
    ), method="wilcox.test", label="p.format") +
    labs(title=title, y="AUC score", x="") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none")
}

p5C <- plot_violin(df, "C3",   "CD14+ mono3 signature")
p5D <- plot_violin(df, "M2",   "M2 Polarization")
p5E <- plot_violin(df, "IL10", "IL-10 signaling")

ggsave(file.path(results_dir,"Fig5C_mono3_sig.png"),  p5C, width=6, height=4)
ggsave(file.path(results_dir,"Fig5D_M2.png"),         p5D, width=6, height=4)
ggsave(file.path(results_dir,"Fig5E_IL10.png"),       p5E, width=6, height=4)

message("✅ Figure 5A–G generated in `results/`")

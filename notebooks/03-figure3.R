# ────────────────────────────────────────────────────────────
# 03-figure3.R: Monocyte sub‐clustering & polarization (Fig 3A–H)
# ────────────────────────────────────────────────────────────

# 1. LIBRARIES
library(Seurat)       # DimPlot, FindMarkers, FetchData
library(dplyr)        # data wrangling
library(ggplot2)      # plotting
library(ggridges)     # ridge density plots
library(cowplot)      # combine panels
library(ggpubr)       # stat_compare_means
library(UCell)        # AddModuleScore_UCell

# 2. I/O
data_file   <- "data/HCCprocesseddata.rds"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# 3. LOAD & SUBSET TO MONOCYTES
se <- readRDS(data_file)
mono_types <- c("CD14+ mono1","CD14+ mono2","CD14+ mono3","CD14+ mono4","CD16+ mono")
Mono <- subset(se, subset = CELLtype %in% mono_types)

# make sure time is ordered and renamed
Mono$time <- factor(recode(Mono$time,"Control"="hc"),
                    levels = c("pre","post","hc"))

# custom colors for subclusters
mono_cols <- c(
  "CD14+ mono1"="#3F007D","CD14+ mono2"="#7A0177",
  "CD14+ mono3"="#C51B8A","CD14+ mono4"="#F768A1",
  "CD16+ mono" ="#FBB4C4"
)

# ────────────────────────────────────────────────────────────
# Fig 3A: UMAP of monocyte subtypes
# ────────────────────────────────────────────────────────────
p3A <- DimPlot(
  Mono, reduction="umap", group.by="CELLtype", cols=mono_cols,
  label=FALSE, pt.size=0.5
) + theme_minimal() +
    theme(legend.text=element_text(size=9),
          axis.text=element_text(size=10))
ggsave(file.path(results_dir,"Fig3A_UMAP_mono.png"),
       p3A, width=5, height=4)

# ────────────────────────────────────────────────────────────
# Fig 3B: (top) density ridges by pseudotime,
#          (bottom) UMAP colored by pseudotime
# (assumes Mono$pseudotime exists in metadata)
# ────────────────────────────────────────────────────────────

meta <- Mono@meta.data %>% select(time, pseudotime, CELLtype)

# ridge densities
p_ridges <- ggplot(meta, aes(x=pseudotime, y=CELLtype, fill=CELLtype)) +
  geom_density_ridges(alpha=0.7) +
  facet_wrap(~time, ncol=3) +
  scale_fill_manual(values=mono_cols) +
  labs(x="Pseudotime", y="") +
  theme_ridges(font_size=10) + theme(legend.position="none")

# split UMAP
p_umap_pt <- DimPlot(
  Mono, reduction="umap", split.by="time", group.by="pseudotime",
  pt.size=0.5
) + scale_color_viridis_c() +
    theme_minimal() + theme(legend.position="none")

# combine
p3B <- plot_grid(p_ridges,p_umap_pt, ncol=1, rel_heights=c(1,1.2))
ggsave(file.path(results_dir,"Fig3B_Pseudotime.png"),
       p3B, width=6, height=8)

# ────────────────────────────────────────────────────────────
# Fig 3C: Box+lines of monocyte‐subtype frequencies over time
# ────────────────────────────────────────────────────────────
dfp <- Mono@meta.data %>%
  count(sample, time, CELLtype) %>%
  group_by(sample, time) %>%
  mutate(pct = n/sum(n)*100) %>%
  ungroup()

p3C <- ggplot(dfp, aes(x=time, y=pct, color=time)) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  geom_line(aes(group=sample), data=filter(dfp,time%in%c("pre","post")),
            color="gray50", alpha=0.5) +
  geom_jitter(width=0.2,size=1) +
  facet_wrap(~CELLtype, nrow=1) +
  stat_compare_means(comparisons=list(c("pre","post")), paired=TRUE) +
  scale_color_manual(values=c("pre"="#F2B342","post"="#5AAA46","hc"="#4F63B0")) +
  labs(y="Percentage (%)", x="") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(results_dir,"Fig3C_CelltypeFreq.png"),
       p3C, width=8, height=3)

# ────────────────────────────────────────────────────────────
# Fig 3D–E: CD14+ mono3 vs CD14+ mono2 volcano & enrichment
# ────────────────────────────────────────────────────────────
# D: volcano
vdat <- FindMarkers(
  Mono, ident.1="CD14+ mono3", ident.2="CD14+ mono2",
  logfc.threshold=0, test.use="MAST", min.pct=0.1
)
vdat$gene <- rownames(vdat)
vdat$p_val_adj[vdat$p_val_adj==0] <- min(vdat$p_val_adj[vdat$p_val_adj>0])*1e-1
vdat$group <- case_when(
  vdat$avg_log2FC>0.5 & vdat$p_val_adj<0.05 ~ "mono3 up",
  vdat$avg_log2FC< -0.5 & vdat$p_val_adj<0.05 ~ "mono2 up",
  TRUE ~ "ns"
)
p3D <- ggplot(vdat, aes(x=avg_log2FC,y=-log10(p_val_adj),color=group)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("mono2 up"=mono_cols["CD14+ mono2"],
                              "mono3 up"=mono_cols["CD14+ mono3"],
                              "ns"="grey70")) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text_repel(data=subset(vdat,group!="ns"),aes(label=gene),
                  size=2, max.overlaps=15) +
  theme_minimal() +
  labs(x="Log2 FC",y="-Log10 FDR",title="CD14+ mono3 vs mono2")
ggsave(file.path(results_dir,"Fig3D_Volcano.png"),
       p3D, width=5, height=4)

# E: bar of top 5 enriched pathways (example hard‐coded)
# replace with your real enrichment table
enr <- data.frame(
  pathway=c("Ubiquitin\nmediated", "Macroautophagy",
            "Antigen\npresentation","Inflammatory\nresponse",
            "Cell–cell\nadhesion"),
  padj = c(2e-9,5e-11,3e-9,2e-6,4e-13)
) %>% mutate(logp=-log10(padj),
             side=ifelse(pathway %in% c("Ubiquitin\nmediated","Macroautophagy"),
                         "down","up"))
p3E <- ggplot(enr,aes(x=logp,y=reorder(pathway,logp),fill=side))+
  geom_col()+
  scale_fill_manual(values=c("down"=mono_cols["CD14+ mono3"],"up"=mono_cols["CD14+ mono2"]))+
  geom_vline(xintercept=0)+
  labs(x="-Log10 padj",y="",title="Pathways enriched")+
  theme_minimal()+theme(legend.position="none")
ggsave(file.path(results_dir,"Fig3E_Enrichment.png"),
       p3E, width=4, height=3)

# ────────────────────────────────────────────────────────────
# Fig 3F–G–H: M1/M2 polarization scores by subtype & time
# ────────────────────────────────────────────────────────────
# define gene lists
M1_genes <- c("IL1B","TNF","NOS2","IL12B","CXCL10")    # trim to your list
M2_genes <- c("ARG1","CD163","MRC1","IL10","TGFB1")    # trim to your list

Mono <- AddModuleScore_UCell(Mono,list(M1_genes), name="M1")
Mono <- AddModuleScore_UCell(Mono,list(M2_genes), name="M2")

md <- FetchData(Mono, vars=c("M1_UCell","M2_UCell","CELLtype","time"))

# F: violin of M1 across subtypes/time
p3F <- ggplot(md, aes(x=time,y=M1_UCell,fill=time))+
  geom_violin(trim=FALSE)+
  stat_compare_means(comparisons=list(c("pre","post"),c("pre","hc")),method="wilcox")+
  facet_wrap(~CELLtype,ncol=5)+
  scale_fill_manual(values=c("pre"="#F2B342","post"="#5AAA46","hc"="#4F63B0"))+
  labs(title="M1 Polarization",y="Score",x="")+
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(results_dir,"Fig3F_M1_violin.png"),
       p3F, width=8, height=4)

# G: violin of M2
p3G <- ggplot(md, aes(x=time,y=M2_UCell,fill=time))+
  geom_violin(trim=FALSE)+
  stat_compare_means(comparisons=list(c("pre","post"),c("pre","hc")),method="wilcox")+
  facet_wrap(~CELLtype,ncol=5)+
  scale_fill_manual(values=c("pre"="#F2B342","post"="#5AAA46","hc"="#4F63B0"))+
  labs(title="M2 Polarization",y="Score",x="")+
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(results_dir,"Fig3G_M2_violin.png"),
       p3G, width=8, height=4)

# H: combined line+box across time for each subtype
p3H <- ggplot(md, aes(x=time,y=M1_UCell,group=CELLtype,color=CELLtype))+
  geom_boxplot(alpha=0.6,outlier.shape=NA)+
  geom_jitter(width=0.15,size=0.8)+
  stat_compare_means(aes(group=time),comparisons=list(c("pre","post"),c("post","hc")),method="wilcox")+
  labs(title="M1 Polarization across subtypes",y="Score",x="")+
  theme_minimal()+theme(legend.position="bottom")
ggsave(file.path(results_dir,"Fig3H_M1_across.png"),
       p3H, width=6, height=4)

message("✅ Figure 3 panels complete. Check results/ for Fig3A–H.")  

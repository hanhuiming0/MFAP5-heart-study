# A ------------------------
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

input.file <- "D:/Dataset/MFAP5/enrichment-kegg-MFAP5-vs-CT1-Up.xls"
code.dir <- "D:/Dataset/MFAP5/code"
figure.dir <- "D:/Dataset/MFAP5/Figures"
date.tag <- "20260914"

dir.create(code.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure.dir, recursive = TRUE, showWarnings = FALSE)

# Read data
dat <- read.delim(
  input.file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    `q-value` = as.numeric(`q-value`),
    `p-value` = as.numeric(`p-value`),
    Enrichment_score = as.numeric(Enrichment_score),
    ListHits = as.numeric(ListHits)
  ) %>%
  filter(
    is.finite(`q-value`),
    is.finite(Enrichment_score)
  )

# Select Top 20 pathways:
# first ranked by q-value (ascending),
# then by Enrichment_score (descending)
plot.dat <- dat %>%
  distinct(Term, .keep_all = TRUE) %>%
  arrange(`q-value`, desc(Enrichment_score)) %>%
  slice_head(n = 20) %>%
  mutate(
    negLog10Q = -log10(pmax(`q-value`, .Machine$double.xmin))
  )

# KEGG category order
category.order <- c(
  "Environmental Information Processing",
  "Human Diseases",
  "Metabolism",
  "Organismal Systems"
)

plot.dat$Classification_level1 <- factor(
  plot.dat$Classification_level1,
  levels = c(
    category.order,
    setdiff(unique(plot.dat$Classification_level1), category.order)
  )
)

# Order pathways within categories
plot.dat <- plot.dat %>%
  arrange(
    Classification_level1,
    Enrichment_score,
    `q-value`
  ) %>%
  mutate(
    Term = factor(Term, levels = unique(Term))
  )

facet.labels <- c(
  "Environmental Information Processing" = "EIP",
  "Human Diseases" = "Human Diseases",
  "Metabolism" = "Metabolism",
  "Organismal Systems" = "Organismal Systems"
)

# Plot
p <- ggplot(
  plot.dat,
  aes(x = Enrichment_score, y = Term)
) +
  geom_point(
    aes(
      size = ListHits,
      color = negLog10Q
    ),
    alpha = 0.95
  ) +
  facet_grid(
    Classification_level1 ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = TRUE,
    labeller = as_labeller(facet.labels)
  ) +
  scale_color_gradientn(
    colors = c(
      "#EFF3FF",
      "#BDD7E7",
      "#6BAED6",
      "#3182BD",
      "#08519C"
    ),
    name = expression(-Log[10] * "(q-value)")
  ) +
  scale_size_continuous(
    range = c(2.7, 7.2),
    name = "Gene count",
    breaks = pretty_breaks(n = 3)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(.03, .08))
  ) +
  labs(
    title = "MFAP5 vs. CT1 (Upregulated genes)",
    subtitle = "KEGG enrichment Top 20",
    x = "Enrichment score",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = .5,
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = .5,
      color = "grey35"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "#E8EDF2",
      linewidth = .4
    ),
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "#E8E8E8",
      color = "#BDBDBD"
    ),
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 8.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 9
    ),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(8, 12, 7, 8)
  )

# Save
base.name <- file.path(
  figure.dir,
  paste0("MFAP5_KEGG_Up_top20_", date.tag)
)

ggsave(
  paste0(base.name, ".pdf"),
  p,
  width = 9.3,
  height = 8.4,
  useDingbats = FALSE
)

write.csv(
  plot.dat %>%
    arrange(`q-value`) %>%
    select(
      id,
      Term,
      Classification_level1,
      Classification_level2,
      ListHits,
      `p-value`,
      `q-value`,
      Enrichment_score
    ),
  paste0(base.name, ".csv"),
  row.names = FALSE
)

# C -----------------
#!/usr/bin/env Rscript

# ============================================================
# Script: MFAP5_Glycolysis_heatmap.R
# Description:
#   Identify significantly altered Glycolysis/Gluconeogenesis
#   genes in MFAP5 vs. CT1 RNA-seq data and generate a heatmap.
#
# Criteria:
#   FDR < 0.05
#   |log2FC| >= log2(1.5)
#
# Input:
#   1. DEG_combinded_nofiltered.xls
#   2. Glycolysis.txt
#
# Output:
#   1. MFAP5_Glycolysis_significant_heatmap.pdf
#   2. MFAP5_Glycolysis_significant_heatmap.png
#   3. MFAP5_Glycolysis_significant_genes_expression.csv
# ============================================================


# -----------------------------
# Load packages
# -----------------------------
suppressPackageStartupMessages({
  library(pheatmap)
})


# -----------------------------
# File paths
# -----------------------------
input.file <- "D:/Dataset/MFAP5/DEG_combinded_nofiltered.xls"
glycolysis.file <- "D:/Dataset/MFAP5/Glycolysis.txt"

code.dir <- "D:/Dataset/MFAP5/codes"
figure.dir <- "D:/Dataset/MFAP5/Figures"

dir.create(
  code.dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure.dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# -----------------------------
# Parameters
# -----------------------------
fc.cutoff <- log2(1.5)
fdr.cutoff <- 0.05

fc.col <- "MFAP5-vs-CT1-all.gene_FC"
p.col <- "MFAP5-vs-CT1-all.gene_pValue"

sample.cols <- c(
  "CT1_1",
  "CT1_2",
  "CT1_3",
  "MFAP5_1",
  "MFAP5_2",
  "MFAP5_3"
)


# -----------------------------
# Read DEG data
# -----------------------------
dat <- read.delim(
  input.file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

stopifnot(
  all(
    c(
      "id",
      fc.col,
      p.col,
      sample.cols
    ) %in% names(dat)
  )
)


# -----------------------------
# Read glycolysis gene list
# -----------------------------
glycolysis <- unique(
  toupper(
    trimws(
      readLines(
        glycolysis.file,
        warn = FALSE
      )
    )
  )
)

glycolysis <- setdiff(
  glycolysis,
  c("", "GLYCOLYSIS")
)


# -----------------------------
# Calculate DEG statistics
# -----------------------------
dat$Gene <- toupper(
  trimws(dat$id)
)

dat$FC <- suppressWarnings(
  as.numeric(dat[[fc.col]])
)

dat$Pvalue <- suppressWarnings(
  as.numeric(dat[[p.col]])
)

dat$FDR <- p.adjust(
  dat$Pvalue,
  method = "BH"
)

dat$log2FC <- ifelse(
  is.finite(dat$FC) & dat$FC > 0,
  log2(dat$FC),
  NA_real_
)


# -----------------------------
# Select significant
# Glycolysis/Gluconeogenesis genes
# -----------------------------
sig <- dat[
  dat$Gene %in% glycolysis &
    is.finite(dat$FDR) &
    dat$FDR < fdr.cutoff &
    abs(dat$log2FC) >= fc.cutoff,
  ,
  drop = FALSE
]

sig <- sig[
  order(
    sig$log2FC,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

sig <- sig[
  !duplicated(sig$Gene),
  ,
  drop = FALSE
]

if (!nrow(sig)) {
  stop(
    "No significant Glycolysis/Gluconeogenesis genes found."
  )
}


# -----------------------------
# Prepare expression matrix
# -----------------------------
expr <- as.matrix(
  data.frame(
    lapply(
      sig[, sample.cols, drop = FALSE],
      as.numeric
    ),
    check.names = FALSE
  )
)

rownames(expr) <- sig$Gene


# Log2 transformation
expr.log <- log2(
  expr + 1
)


# Row-wise Z-score normalization
expr.z <- t(
  scale(
    t(expr.log)
  )
)

expr.z[
  !is.finite(expr.z)
] <- 0


# Limit Z-score range
expr.z <- pmax(
  pmin(
    expr.z,
    2
  ),
  -2
)


# -----------------------------
# Column annotation
# -----------------------------
annotation.col <- data.frame(
  Group = factor(
    c(
      rep("CT1", 3),
      rep("MFAP5", 3)
    ),
    levels = c(
      "CT1",
      "MFAP5"
    )
  ),
  row.names = sample.cols
)


# -----------------------------
# Row annotation
# -----------------------------
annotation.row <- data.frame(
  Regulation = factor(
    ifelse(
      sig$log2FC > 0,
      "Up",
      "Down"
    ),
    levels = c(
      "Up",
      "Down"
    )
  ),
  row.names = sig$Gene
)


# -----------------------------
# Annotation colors
# -----------------------------
annotation.colors <- list(
  Group = c(
    CT1 = "#17becf",
    MFAP5 = "#ff9896"
  ),
  Regulation = c(
    Up = "#D95F5F",
    Down = "#4C78A8"
  )
)


# -----------------------------
# Heatmap colors
# -----------------------------
heat.colors <- colorRampPalette(
  c(
    "#2B4C9B",
    "#F7F7F7",
    "#D95F6A"
  )
)(101)


# -----------------------------
# Heatmap function
# -----------------------------
plot.heatmap <- function(
  filename,
  type
) {

  if (type == "pdf") {

    pdf(
      filename,
      width = 6.2,
      height = max(
        8,
        0.25 * nrow(expr.z) + 2.3
      ),
      useDingbats = FALSE
    )

  } else {

    png(
      filename,
      width = 2200,
      height = max(
        2800,
        110 * nrow(expr.z) + 700
      ),
      res = 320
    )
  }

  pheatmap(
    expr.z,
    color = heat.colors,
    breaks = seq(
      -2,
      2,
      length.out = 102
    ),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation.col,
    annotation_row = annotation.row,
    annotation_colors = annotation.colors,
    annotation_names_row = FALSE,
    gaps_col = 3,
    border_color = "white",
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = 11,
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = 0,
    treeheight_col = 0,
    main = paste0(
      "Significant Glycolysis / Gluconeogenesis genes\n",
      "MFAP5 vs. CT1"
    ),
    legend_breaks = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    legend_labels = c(
      "-2",
      "-1",
      "0",
      "1",
      "2"
    )
  )

  dev.off()
}


# -----------------------------
# Export heatmaps
# -----------------------------
plot.heatmap(
  file.path(
    figure.dir,
    "MFAP5_Glycolysis_significant_heatmap.pdf"
  ),
  "pdf"
)

plot.heatmap(
  file.path(
    figure.dir,
    "MFAP5_Glycolysis_significant_heatmap.png"
  ),
  "png"
)


# -----------------------------
# Export significant genes
# -----------------------------
out <- data.frame(
  Gene = sig$Gene,
  FC = sig$FC,
  log2FC = sig$log2FC,
  Pvalue = sig$Pvalue,
  FDR = sig$FDR,
  Regulation = ifelse(
    sig$log2FC > 0,
    "Up",
    "Down"
  ),
  expr,
  check.names = FALSE
)

write.csv(
  out,
  file.path(
    figure.dir,
    "MFAP5_Glycolysis_significant_genes_expression.csv"
  ),
  row.names = FALSE
)


# -----------------------------
# Summary
# -----------------------------
message(
  "Completed: ",
  nrow(sig),
  " significant Glycolysis/Gluconeogenesis genes (",
  sum(sig$log2FC > 0),
  " up, ",
  sum(sig$log2FC < 0),
  " down)."
)

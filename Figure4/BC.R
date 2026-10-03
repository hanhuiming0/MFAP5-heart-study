# B ---------------------------

#!/usr/bin/env Rscript

# ============================================================
# Script: MFAP5_ITGAV_spatial_colocalization.R
#
# Description:
#   Evaluate the spatial co-localization of MFAP5 and ITGAV
#   in human cardiac spatial transcriptomic datasets.
#
#   For each gene, the 20 spots with the highest expression
#   are defined as high-expression spots. Spatial spots are
#   subsequently classified into four groups:
#
#     1. Both low
#     2. MFAP5 high / ITGAV low
#     3. MFAP5 low / ITGAV high
#     4. Both high
#
#   RCTD-derived cardiomyocyte proportions are additionally
#   visualized for comparison with gene co-localization.
#
# Output:
#   One PDF per spatial transcriptomic sample containing:
#     - Cardiomyocyte proportion
#     - MFAP5 expression
#     - ITGAV expression
#     - MFAP5/ITGAV co-localization
# ============================================================


# -----------------------------
# Load packages
# -----------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})


# -----------------------------
# Parameters
# -----------------------------
pair <- c(
  "MFAP5",
  "ITGAV"
)

top.n <- 20

pt.size <- 0.1


# -----------------------------
# File paths
# -----------------------------
data.dir <- paste0(
  "/storage/data/Hhm/Heart_Study/",
  "data/ST/3.heartcellatlas/Data_LV.SP_r"
)

rctd.dir <- paste0(
  "/storage/data/Hhm/Heart_Study/NZH"
)

output.dir <- rctd.dir

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# Sample information
# -----------------------------
sample_info <- data.frame(
  sample_id = c(
    "HCAHeartST11350375",
    "HCAHeartST10550732",
    "HCAHeartST8795938",
    "HCAHeartST8795939"
  ),

  seurat_file = c(
    "SP_AH1_HCAHeartST11350375.h5ad_seurat_obj.rds",
    "SP_D3_HCAHeartST10550732.h5ad_seurat_obj.rds",
    "SP_D5_HCAHeartST8795938.h5ad_seurat_obj.rds",
    "SP_D5_HCAHeartST8795939.h5ad_seurat_obj.rds"
  ),

  rctd_file = c(
    "SP_AH1_HCAHeartST11350375.h5ad_RCTD_Result.RDS",
    "SP_D3_HCAHeartST10550732.h5ad_RCTD_Result.RDS",
    "SP_D5_HCAHeartST8795938.h5ad_RCTD_Result.RDS",
    "SP_D5_HCAHeartST8795939.h5ad_RCTD_Result.RDS"
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# Process each spatial transcriptomic sample
# ============================================================

for (i in seq_len(nrow(sample_info))) {

  # ---------------------------
  # Sample information
  # ---------------------------
  sample <- sample_info[i, "sample_id"]

  message(
    "Processing: ",
    sample
  )


  # ---------------------------
  # Load Seurat object
  # ---------------------------
  ST.data <- readRDS(
    file.path(
      data.dir,
      sample_info[i, "seurat_file"]
    )
  )


  # ---------------------------
  # Normalize expression
  # ---------------------------
  ST.data <- NormalizeData(
    ST.data,
    verbose = FALSE
  )


  # ---------------------------
  # Spatial coordinates
  # ---------------------------
  spatial_coord <- data.frame(
    GetTissueCoordinates(
      ST.data
    )
  )

  spatial_coord$barcodeID <-
    rownames(spatial_coord)

  location <- spatial_coord[
    ,
    c(
      "imagerow",
      "imagecol"
    )
  ]


  # ---------------------------
  # Extract MFAP5/ITGAV expression
  # ---------------------------
  expr <- FetchData(
    ST.data,
    vars = pair
  )

  ncell <- nrow(expr)


  # ---------------------------
  # Load RCTD results
  # ---------------------------
  ST.RCTD <- readRDS(
    file.path(
      rctd.dir,
      sample_info[i, "rctd_file"]
    )
  )


  # ---------------------------
  # Normalize RCTD weights
  # ---------------------------
  weights <- ST.RCTD@results$weights

  norm_weights <- normalize_weights(
    weights
  )

  norm_weights <- as.data.frame(
    norm_weights
  )


  # Add cardiomyocyte proportion
  ST.data$Cardiomyocyte <-
    norm_weights$Cardiomyocyte[
      match(
        colnames(ST.data),
        rownames(norm_weights)
      )
    ]


  # ==========================================================
  # Spatial feature plots
  # ==========================================================

  p1 <- SpatialFeaturePlot(
    ST.data,
    features = "Cardiomyocyte"
  )


  p2 <- SpatialFeaturePlot(
    ST.data,
    features = pair[1]
  )


  p3 <- SpatialFeaturePlot(
    ST.data,
    features = pair[2]
  )


  # ==========================================================
  # Define high-expression spots
  # ==========================================================

  gene1 <- expr[, pair[1]]
  gene2 <- expr[, pair[2]]


  # Number of spots to select
  n1.select <- min(
    top.n,
    sum(gene1 > 0)
  )

  n2.select <- min(
    top.n,
    sum(gene2 > 0)
  )


  # Select Top 20 expressing spots
  n1 <- order(
    gene1,
    decreasing = TRUE
  )

  n1 <- n1[
    gene1[n1] > 0
  ]

  n1 <- head(
    n1,
    n1.select
  )


  n2 <- order(
    gene2,
    decreasing = TRUE
  )

  n2 <- n2[
    gene2[n2] > 0
  ]

  n2 <- head(
    n2,
    n2.select
  )


  # ==========================================================
  # Classify spots into four groups
  # ==========================================================

  expcol <- rep(
    0,
    ncell
  )


  # MFAP5 high
  expcol[n1] <- 1


  # ITGAV high
  expcol[n2] <- 2


  # Both high
  expcol[
    intersect(
      n1,
      n2
    )
  ] <- 3


  # ---------------------------
  # Prepare plotting data
  # ---------------------------
  tmp <- data.frame(
    x = location[, 1],
    y = location[, 2],

    Exp = factor(
      expcol,
      levels = 0:3,
      labels = c(
        "Both low",
        "MFAP5 high",
        "ITGAV high",
        "Both high"
      )
    )
  )


  # ==========================================================
  # Co-localization plot
  # ==========================================================

  p4 <- ggplot(
    tmp,
    aes(
      x = x,
      y = y,
      color = Exp
    )
  ) +

    geom_point(
      size = pt.size
    ) +

    scale_color_manual(
      values = c(
        "Both low" = "#324E97",
        "MFAP5 high" = "#F6C46A",
        "ITGAV high" = "#396C1E",
        "Both high" = "#A22B21"
      ),

      guide = guide_legend(
        override.aes = list(
          size = 4
        )
      )
    ) +

    ggtitle(
      paste(
        pair,
        collapse = "_"
      )
    ) +

    scale_x_reverse() +

    coord_flip() +

    theme_minimal() +

    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )


  # ==========================================================
  # Combine plots
  # ==========================================================

  combined.plot <-
    p1 +
    p2 +
    p3 +
    p4 +
    plot_layout(
      ncol = 4
    )


  # ==========================================================
  # Save PDF
  # ==========================================================

  output.file <- file.path(
    output.dir,
    paste0(
      sample,
      ".MFAP5_ITGAV.colocation.pdf"
    )
  )


  pdf(
    output.file,
    width = 8.6,
    height = 3,
    useDingbats = FALSE
  )

  print(
    combined.plot
  )

  dev.off()


  message(
    "  MFAP5 high spots: ",
    length(n1),
    "; ITGAV high spots: ",
    length(n2),
    "; Both high spots: ",
    length(
      intersect(
        n1,
        n2
      )
    )
  )
}


message(
  "Spatial co-localization analysis completed."
)












# C ----------------------------------------
#!/usr/bin/env Rscript

# ============================================================
# Script: GSE249925_MFAP5_ITGAV_correlation.R
#
# Description:
#   Evaluate the correlation between MFAP5 and ITGAV
#   expression in HCM samples from GSE249925.
#
# Statistical analysis:
#   Pearson correlation
#   Expression transformation: log10(x + 1)
#
# Output:
#   GSE249925_MFAP5_ITGAV_correlation.pdf
# ============================================================


# -----------------------------
# Load packages
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(psych)
})


# -----------------------------
# File paths
# -----------------------------
expr.file <- paste0(
  "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
  "GSE249925/Result/expr_symbol2.RDS"
)

deg.file <- paste0(
  "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
  "GSE249925/Result/GSE249925_DEG.RDS"
)

sample.file <- paste0(
  "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
  "GSE249925/GSE249925_Sample_Info.txt"
)

output.file <- paste0(
  "E:/1-课题/其他/MFAP5/表达相关性/",
  "GSE249925_MFAP5-ITGAV-表达相关性.pdf"
)


# -----------------------------
# Read expression and DEG data
# -----------------------------
GSE249925.expr_symbol <- readRDS(
  expr.file
)

GSE249925.DEGs <- readRDS(
  deg.file
)


# Optional:
# inspect MFAP5 differential-expression result
print(
  GSE249925.DEGs[
    "MFAP5",
    ,
    drop = FALSE
  ]
)


# -----------------------------
# Read sample information
# -----------------------------
GSE249925.Sample_Info <- read.delim(
  sample.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Rename phenotype column
colnames(
  GSE249925.Sample_Info
)[3] <- "Sample_Type"


# Rename Control as Normal
GSE249925.Sample_Info$Sample_Type[
  GSE249925.Sample_Info$Sample_Type == "Control"
] <- "Normal"


# -----------------------------
# Match sample IDs
# -----------------------------
rownames(
  GSE249925.Sample_Info
) <- GSE249925.Sample_Info$Sample_title


# Match R-generated column names beginning with numbers
rownames(
  GSE249925.Sample_Info
)[1:91] <- paste0(
  "X",
  GSE249925.Sample_Info$Sample_title[1:91]
)


# -----------------------------
# Retain Normal and HCM samples
# -----------------------------
GSE249925.Sample_Info <- GSE249925.Sample_Info[
  GSE249925.Sample_Info$Sample_Type %in%
    c(
      "Normal",
      "HCM"
    ),
  ,
  drop = FALSE
]


# Subset expression matrix
GSE249925.expr_symbol <- GSE249925.expr_symbol[
  ,
  rownames(
    GSE249925.Sample_Info
  ),
  drop = FALSE
]


# -----------------------------
# Check target genes
# -----------------------------
target.genes <- c(
  "MFAP5",
  "ITGAV"
)

if (
  !all(
    target.genes %in%
      rownames(
        GSE249925.expr_symbol
      )
  )
) {
  stop(
    "MFAP5 or ITGAV was not found in the expression matrix."
  )
}


# -----------------------------
# Prepare expression data
# -----------------------------
plot.data <- data.frame(
  Sample_ID = rownames(
    GSE249925.Sample_Info
  ),

  Sample_Type =
    GSE249925.Sample_Info$Sample_Type,

  MFAP5 = as.numeric(
    GSE249925.expr_symbol[
      "MFAP5",
      rownames(
        GSE249925.Sample_Info
      )
    ]
  ),

  ITGAV = as.numeric(
    GSE249925.expr_symbol[
      "ITGAV",
      rownames(
        GSE249925.Sample_Info
      )
    ]
  )
)


# -----------------------------
# Log transformation
# -----------------------------
plot.data$MFAP5 <- log10(
  plot.data$MFAP5 + 1
)

plot.data$ITGAV <- log10(
  plot.data$ITGAV + 1
)


# -----------------------------
# Retain HCM samples
# -----------------------------
plot.data <- plot.data[
  plot.data$Sample_Type == "HCM",
  ,
  drop = FALSE
]


# ============================================================
# Pearson correlation
# ============================================================

myCor <- corr.test(
  plot.data$MFAP5,
  plot.data$ITGAV,
  use = "pairwise",
  method = "pearson",
  adjust = "none"
)


# Correlation label
label_text <- paste0(
  "R = ",
  round(
    myCor$r,
    2
  ),
  ", P = ",
  format(
    myCor$p,
    scientific = TRUE,
    digits = 2
  )
)


# ============================================================
# Correlation plot
# ============================================================

p <- ggplot(
  plot.data,
  aes(
    x = MFAP5,
    y = ITGAV
  )
) +

  # Individual HCM samples
  geom_point(
    size = 4,
    fill = "#ff8822",
    color = "black",
    shape = 21,
    stroke = 1
  ) +

  # Linear regression and 95% confidence interval
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "#F29325",
    fill = "#aaaaaa",
    linewidth = 0.75
  ) +

  # Pearson correlation
  annotate(
    "text",
    x = 2.9,
    y = 4.5,
    label = label_text,
    size = 6
  ) +

  labs(
    x = expression(
      MFAP5~expression
    ),
    y = expression(
      ITGAV~expression
    )
  ) +

  theme_bw() +

  theme(
    panel.border =
      element_blank(),

    panel.grid =
      element_line(
        color = "#e0e0e0"
      ),

    axis.title =
      element_text(
        size = 18
      ),

    axis.text =
      element_text(
        size = 16
      ),

    axis.ticks =
      element_blank(),

    plot.margin =
      margin(
        10,
        10,
        10,
        10
      )
  )


# -----------------------------
# Save figure
# -----------------------------
ggsave(
  output.file,
  p,
  width = 4,
  height = 4,
  useDingbats = FALSE
)


# -----------------------------
# Print statistics
# -----------------------------
message(
  "Pearson correlation completed."
)

message(
  "R = ",
  round(
    myCor$r,
    3
  ),
  "; P = ",
  signif(
    myCor$p,
    3
  ),
  "; n = ",
  nrow(
    plot.data
  )
)

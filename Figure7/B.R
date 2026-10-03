# B------------------
#!/usr/bin/env Rscript

# ============================================================
# Script: cross_dataset_TF_analysis.R
#
# Description:
#   Cross-dataset differential expression analysis of selected
#   transcription factors across human HCM/HF datasets and
#   mouse TAC datasets.
#
# Analysis strategy:
#   - Raw count matrices:
#       edgeR TMM normalization + limma-voom
#
#   - Continuous expression matrices:
#       limma after dataset-specific transformation
#
#   - Differential expression direction:
#       Disease/TAC minus Control/Sham
#
#   - Output:
#       1. Per-dataset differential-expression tables
#       2. Per-dataset TF results
#       3. Dataset-selection audit
#       4. Analysis manifest
#       5. Integrated TF results
#       6. Cross-dataset bubble heatmap
#       7. Analysis README
#
# Usage:
#   Rscript cross_dataset_TF_analysis.R <output_directory>
#
# Example:
#   Rscript cross_dataset_TF_analysis.R results/TF_analysis
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})


# ============================================================
# 2. Output directory
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)

outdir <- args[1]

if (
  is.na(outdir) ||
  !nzchar(outdir)
) {
  stop(
    paste0(
      "Output directory required.\n",
      "Usage: Rscript cross_dataset_TF_analysis.R <output_directory>"
    )
  )
}


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    outdir,
    "per_dataset"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Candidate transcription factors
# ============================================================

tfs <- c(
  "PRDM4",
  "NR3C1",
  "KLF12",
  "STAT3",
  "MAFG",
  "KLF4",
  "NR4A2",
  "CREB1",
  "HNF1B",
  "EBF1",
  "SP1",
  "TBX6",
  "SPI1",
  "NFE2L2",
  "BHLHA15",
  "KLF1",
  "NFATC2",
  "DLX5",
  "ATOH1",
  "IRF1",
  "VDR"
)


manifest <- list()
integrated <- list()


# ============================================================
# 4. Helper functions
# ============================================================


# ------------------------------------------------------------
# Collapse rows sharing the same gene symbol
# ------------------------------------------------------------

collapse_symbols <- function(
  mat,
  symbols
) {

  keep <- !is.na(symbols) &
    nzchar(symbols)

  mat <- as.matrix(
    mat[
      keep,
      ,
      drop = FALSE
    ]
  )

  symbols <- symbols[
    keep
  ]

  rowsum(
    mat,
    group = symbols,
    reorder = FALSE
  )
}


# ------------------------------------------------------------
# Map Entrez IDs to gene symbols
# ------------------------------------------------------------

map_entrez <- function(
  ids,
  species
) {

  db <- if (
    species == "human"
  ) {
    org.Hs.eg.db
  } else {
    org.Mm.eg.db
  }

  unname(
    mapIds(
      db,
      keys = as.character(ids),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  )
}


# ------------------------------------------------------------
# Finalize and export one differential-expression result
# ------------------------------------------------------------

finish <- function(
  tab,
  dataset,
  species,
  comparison,
  source,
  n,
  nd,
  nc,
  method
) {

  tab$dataset <- dataset
  tab$species <- species
  tab$comparison <- comparison

  tab <- tab[
    ,
    c(
      "dataset",
      "species",
      "comparison",
      "gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]


  # Export complete differential-expression table
  write.csv(
    tab,
    file.path(
      outdir,
      "per_dataset",
      paste0(
        dataset,
        "_all_DE.csv"
      )
    ),
    row.names = FALSE
  )


  # Select candidate TFs
  want <- toupper(
    tab$gene
  ) %in% tfs

  tf <- tab[
    want,
    ,
    drop = FALSE
  ]

  tf$gene <- toupper(
    tf$gene
  )


  # Preserve predefined TF order
  tf <- tf[
    match(
      tfs,
      tf$gene
    ),
    ,
    drop = FALSE
  ]


  # Fill missing TFs
  missing <- is.na(
    tf$gene
  )

  tf$gene[
    missing
  ] <- tfs[
    missing
  ]

  tf$dataset[
    missing
  ] <- dataset

  tf$species[
    missing
  ] <- species

  tf$comparison[
    missing
  ] <- comparison


  # Nominal P-value significance
  tf$significance <- ifelse(
    is.na(tf$P.Value),
    "NA",
    ifelse(
      tf$P.Value < 0.001,
      "***",
      ifelse(
        tf$P.Value < 0.01,
        "**",
        ifelse(
          tf$P.Value < 0.05,
          "*",
          "ns"
        )
      )
    )
  )


  # Export TF-specific result
  write.csv(
    tf,
    file.path(
      outdir,
      "per_dataset",
      paste0(
        dataset,
        "_TF_results.csv"
      )
    ),
    row.names = FALSE
  )


  integrated[
    [dataset]
  ] <<- tf


  manifest[
    [dataset]
  ] <<- data.frame(
    dataset = dataset,
    species = species,
    comparison = comparison,
    n_total = n,
    n_disease = nd,
    n_control = nc,
    genes_tested = nrow(tab),
    method = method,
    source_file = source,
    status = "Included",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Differential expression for count matrices
# ------------------------------------------------------------

run_counts <- function(
  mat,
  group,
  dataset,
  species,
  comparison,
  source
) {

  mat <- round(
    as.matrix(mat)
  )

  storage.mode(
    mat
  ) <- "numeric"


  stopifnot(
    ncol(mat) == length(group),
    all(
      group %in%
        c(
          "Control",
          "Disease"
        )
    )
  )


  # CPM filtering
  keep <- rowSums(
    cpm(mat) >= 1
  ) >= min(
    table(group)
  )


  y <- DGEList(
    mat[
      keep,
      ,
      drop = FALSE
    ],
    group = factor(
      group,
      levels = c(
        "Control",
        "Disease"
      )
    )
  )


  # TMM normalization
  y <- calcNormFactors(
    y
  )


  # Design matrix
  design <- model.matrix(
    ~factor(
      group,
      levels = c(
        "Control",
        "Disease"
      )
    )
  )


  # voom + limma
  v <- voom(
    y,
    design,
    plot = FALSE
  )


  fit <- lmFit(
    v,
    design
  )

  fit <- eBayes(
    fit,
    robust = TRUE
  )


  tab <- topTable(
    fit,
    coef = 2,
    number = Inf,
    sort.by = "none"
  )


  tab$gene <- if (
    "ID" %in%
      names(tab)
  ) {
    as.character(
      tab$ID
    )
  } else {
    rownames(tab)
  }


  rownames(
    tab
  ) <- NULL


  tab <- tab[
    ,
    c(
      "gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]


  finish(
    tab = tab,
    dataset = dataset,
    species = species,
    comparison = comparison,
    source = source,
    n = ncol(mat),
    nd = sum(
      group == "Disease"
    ),
    nc = sum(
      group == "Control"
    ),
    method = "counts:TMM+voom-limma"
  )
}


# ------------------------------------------------------------
# Differential expression for continuous matrices
# ------------------------------------------------------------

run_continuous <- function(
  mat,
  group,
  dataset,
  species,
  comparison,
  source,
  transform = "none"
) {

  mat <- as.matrix(
    mat
  )

  storage.mode(
    mat
  ) <- "numeric"


  if (
    transform == "log2p1"
  ) {
    mat <- log2(
      mat + 1
    )
  }


  if (
    transform == "log2p01"
  ) {
    mat <- log2(
      mat + 0.1
    )
  }


  stopifnot(
    ncol(mat) == length(group),
    all(
      group %in%
        c(
          "Control",
          "Disease"
        )
    )
  )


  design <- model.matrix(
    ~factor(
      group,
      levels = c(
        "Control",
        "Disease"
      )
    )
  )


  fit <- lmFit(
    mat,
    design
  )

  fit <- eBayes(
    fit,
    robust = TRUE
  )


  tab <- topTable(
    fit,
    coef = 2,
    number = Inf,
    sort.by = "none"
  )


  tab$gene <- if (
    "ID" %in%
      names(tab)
  ) {
    as.character(
      tab$ID
    )
  } else {
    rownames(tab)
  }


  rownames(
    tab
  ) <- NULL


  tab <- tab[
    ,
    c(
      "gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]


  finish(
    tab = tab,
    dataset = dataset,
    species = species,
    comparison = comparison,
    source = source,
    n = ncol(mat),
    nd = sum(
      group == "Disease"
    ),
    nc = sum(
      group == "Control"
    ),
    method = paste0(
      "continuous:",
      transform,
      "+limma"
    )
  )
}


# ============================================================
# 5. Human datasets
# ============================================================


# ------------------------------------------------------------
# GSE36961
# Human HCM vs control
# ------------------------------------------------------------

f <- paste0(
  "D:/Dataset/2026 Heart/Hs/Bulk/HCM/",
  "GSE36961/GSE36961_expression.rds"
)

x <- readRDS(
  f
)

info <- read.delim(
  paste0(
    "D:/Dataset/2026 Heart/Hs/Bulk/HCM/",
    "GSE36961/GSE36961_Sample_Info.txt"
  ),
  check.names = FALSE
)


vcols <- grep(
  "_mRNA$",
  names(x),
  value = TRUE
)

keys <- sub(
  "_mRNA$",
  "",
  vcols
)

title_key <- sub(
  " \\(mRNA\\)$",
  "",
  info$Sample_title
)

idx <- match(
  keys,
  title_key
)


keep <- !is.na(idx) &
  info$Disease_State[
    idx
  ] %in%
  c(
    "control",
    "hypertrophic cardiomyopathy (HCM)"
  )


mat <- as.matrix(
  x[
    ,
    vcols[
      keep
    ],
    drop = FALSE
  ]
)

rownames(
  mat
) <- x$ID_REF


group <- ifelse(
  info$Disease_State[
    idx[
      keep
    ]
  ] == "control",
  "Control",
  "Disease"
)


run_continuous(
  mat = mat,
  group = group,
  dataset = "GSE36961_HCM",
  species = "Homo sapiens",
  comparison = "HCM vs control",
  source = f,
  transform = "log2p1"
)


# ------------------------------------------------------------
# GSE133054
# Human HCM and HF vs NCM
# ------------------------------------------------------------

f <- paste0(
  "D:/Dataset/2026 Heart/Hs/Bulk/HCM-HF/",
  "GSE133054/RAW/GSE133054_human_heart_tissue_readCount.csv.gz"
)

x <- read.csv(
  gzfile(f),
  check.names = FALSE
)


genes <- map_entrez(
  x[
    [1]
  ],
  "human"
)

mat <- collapse_symbols(
  x[
    ,
    -1
  ],
  genes
)


for (
  disease in c(
    "HCM",
    "HF"
  )
) {

  cols <- c(
    grep(
      paste0(
        "^",
        disease
      ),
      colnames(mat)
    ),
    grep(
      "^NCM",
      colnames(mat)
    )
  )


  grp <- ifelse(
    grepl(
      "^NCM",
      colnames(mat)[
        cols
      ]
    ),
    "Control",
    "Disease"
  )


  run_counts(
    mat = mat[
      ,
      cols,
      drop = FALSE
    ],
    group = grp,
    dataset = paste0(
      "GSE133054_",
      disease
    ),
    species = "Homo sapiens",
    comparison = paste0(
      disease,
      " vs normal"
    ),
    source = f
  )
}


# ------------------------------------------------------------
# GSE141910
# Human HCM vs Normal
# log2 TPM
# ------------------------------------------------------------

f <- paste0(
  "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
  "GSE141910/GSE141910.log2_tpm_matrix.rds"
)

mat <- readRDS(
  f
)

info <- read.delim(
  paste0(
    "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
    "GSE141910/GSE141910_Sample_Info_add2.txt"
  ),
  check.names = FALSE
)


gsm <- sub(
  "_.*$",
  "",
  colnames(mat)
)

idx <- match(
  gsm,
  info$Sample_geo_accession
)


keep <- !is.na(idx) &
  info$Sample_Type[
    idx
  ] %in%
  c(
    "Normal",
    "HCM"
  )


run_continuous(
  mat = mat[
    ,
    keep,
    drop = FALSE
  ],
  group = ifelse(
    info$Sample_Type[
      idx[
        keep
      ]
    ] == "Normal",
    "Control",
    "Disease"
  ),
  dataset = "GSE141910_HCM",
  species = "Homo sapiens",
  comparison = "HCM vs normal",
  source = f,
  transform = "none"
)


# ------------------------------------------------------------
# Human count datasets with sample metadata
# ------------------------------------------------------------

human_count_sets <- list(

  GSE160997_HCM = list(
    expr = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE160997/Result/expr_symbol2.RDS"
    ),
    info = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE160997/GSE160997_Sample_Info.txt"
    ),
    sample = "Sample_geo_accession",
    group = "Sample_Type",
    control = "Normal",
    disease = "HCM"
  ),

  GSE180313_HCM = list(
    expr = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE180313/Result/expr_symbol2.RDS"
    ),
    info = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE180313/GSE180313_Sample_Info.txt"
    ),
    sample = "Sample_title",
    group = "Sample_Type",
    control = "Donor",
    disease = "HCM"
  ),

  GSE249925_HCM = list(
    expr = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE249925/Result/expr_symbol2.RDS"
    ),
    info = paste0(
      "D:/Dataset/HF/Homo sapiens/Bulk/Train/",
      "GSE249925/GSE249925_Sample_Info.txt"
    ),
    sample = "Sample_title",
    group = "Group",
    control = "Control",
    disease = "HCM"
  )
)


for (
  nm in names(
    human_count_sets
  )
) {

  z <- human_count_sets[
    [nm]
  ]

  mat <- readRDS(
    z$expr
  )

  info <- read.delim(
    z$info,
    check.names = FALSE
  )


  sample_keys <- if (
    nm == "GSE249925_HCM"
  ) {
    make.names(
      info[
        [z$sample]
      ]
    )
  } else {
    info[
      [z$sample]
    ]
  }


  idx <- match(
    colnames(mat),
    sample_keys
  )


  keep <- !is.na(idx) &
    info[
      [z$group]
    ][
      idx
    ] %in%
    c(
      z$control,
      z$disease
    )


  run_counts(
    mat = mat[
      ,
      keep,
      drop = FALSE
    ],
    group = ifelse(
      info[
        [z$group]
      ][
        idx[
          keep
        ]
      ] == z$control,
      "Control",
      "Disease"
    ),
    dataset = nm,
    species = "Homo sapiens",
    comparison = "HCM vs normal",
    source = z$expr
  )
}


# ------------------------------------------------------------
# GSE89714
# Human hypertrophic vs non-hypertrophic heart
# RPKM
# ------------------------------------------------------------

f <- paste0(
  "D:/Dataset/Heart Failure/Bulk/",
  "GSE89714/GSE89714_Annotation.txt"
)

x <- read.delim(
  f,
  check.names = FALSE,
  na.strings = c(
    "-",
    "NA"
  )
)


genes <- x[
  [1]
]

mat <- as.matrix(
  x[
    ,
    -1
  ]
)

mat[
  is.na(mat)
] <- 0

rownames(
  mat
) <- genes


run_continuous(
  mat = mat,
  group = ifelse(
    grepl(
      "^HY",
      colnames(mat)
    ),
    "Disease",
    "Control"
  ),
  dataset = "GSE89714_HCM",
  species = "Homo sapiens",
  comparison = "hypertrophic vs non-hypertrophic heart",
  source = f,
  transform = "log2p01"
)


# ============================================================
# 6. Mouse TAC datasets
# ============================================================


# ------------------------------------------------------------
# GSE133054
# Mouse TAC time-course pooled vs sham
# ------------------------------------------------------------

f <- paste0(
  "D:/Dataset/2026 Heart/Hs/Bulk/HCM-HF/",
  "GSE133054/RAW/GSE133054_mouse_TAC_bulk_CM_readCount.csv.gz"
)

x <- read.csv(
  gzfile(f),
  check.names = FALSE
)


genes <- map_entrez(
  x[
    [1]
  ],
  "mouse"
)

mat <- collapse_symbols(
  x[
    ,
    -1
  ],
  genes
)


run_counts(
  mat = mat,
  group = ifelse(
    grepl(
      "sham",
      colnames(mat),
      ignore.case = TRUE
    ),
    "Control",
    "Disease"
  ),
  dataset = "GSE133054_TAC",
  species = "Mus musculus",
  comparison = "TAC (2/5/8/11 weeks pooled) vs sham",
  source = f
)


# ------------------------------------------------------------
# Other mouse TAC datasets
# ------------------------------------------------------------

mouse_sets <- list(

  GSE228199_TAC = list(
    f = paste0(
      "D:/Dataset/HF/Mus musculus/Bulk/",
      "GSE228199/GSE228199.count.RDS"
    ),
    grp = c(
      rep(
        "Disease",
        3
      ),
      rep(
        "Control",
        3
      )
    ),
    cmp = "TAC 6 weeks vs sham"
  ),

  GSE235601_TAC = list(
    f = paste0(
      "D:/Dataset/HF/Mus musculus/Bulk/",
      "GSE235601/GSE235601.data.RDS"
    ),
    grp = c(
      rep(
        "Control",
        4
      ),
      rep(
        "Disease",
        4
      )
    ),
    cmp = "WT TAC 5 weeks vs WT sham"
  ),

  GSE247309_TAC = list(
    f = paste0(
      "D:/Dataset/HF/Mus musculus/Bulk/",
      "GSE247309/GSE247309.data.RDS"
    ),
    grp = NULL,
    cmp = "TAC saline 3 weeks vs sham saline"
  ),

  GSE262105_TAC = list(
    f = paste0(
      "D:/Dataset/HF/Mus musculus/Bulk/",
      "GSE262105/GSE262105.data.RDS"
    ),
    grp = c(
      rep(
        "Disease",
        3
      ),
      rep(
        "Control",
        6
      )
    ),
    cmp = "TAC 8 weeks vs sham"
  ),

  GSE262894_TAC = list(
    f = paste0(
      "D:/Dataset/HF/Mus musculus/Bulk/",
      "GSE262894/GSE262894.data.RDS"
    ),
    grp = c(
      rep(
        "Control",
        3
      ),
      rep(
        "Disease",
        3
      )
    ),
    cmp = "control-genotype TAC 8 weeks vs sham"
  )
)


for (
  nm in names(
    mouse_sets
  )
) {

  z <- mouse_sets[
    [nm]
  ]

  mat <- readRDS(
    z$f
  )

  grp <- z$grp


  if (
    is.null(grp)
  ) {

    grp <- ifelse(
      grepl(
        "_S(2|12|16)$",
        colnames(mat)
      ),
      "Disease",
      "Control"
    )
  }


  run_counts(
    mat = mat,
    group = grp,
    dataset = nm,
    species = "Mus musculus",
    comparison = z$cmp,
    source = z$f
  )
}


# ============================================================
# 7. Dataset inclusion/exclusion audit
# ============================================================

selection <- data.frame(

  dataset = c(
    names(manifest),
    "GSE130036",
    "GSE36074",
    "GSE56348",
    "GSE245034",
    "GSE249409",
    "GSE249411",
    "GSE255532"
  ),

  decision = c(
    rep(
      "Included",
      length(manifest)
    ),
    rep(
      "Excluded",
      7
    )
  ),

  reason = c(

    rep(
      "Eligible contrast and complete local files",
      length(manifest)
    ),

    "Existing QC note: too few DE genes; excluded",

    "User explicitly excluded",

    "User explicitly excluded",

    "HF label but no explicit TAC contrast in metadata summary",

    "HFD+mTAC is confounded, not plain TAC",

    "HF label but no explicit TAC contrast in metadata summary",

    "ISO model, not TAC"
  ),

  stringsAsFactors = FALSE
)


write.csv(
  selection,
  file.path(
    outdir,
    "00_dataset_selection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. Analysis manifest
# ============================================================

mani <- do.call(
  rbind,
  manifest
)


write.csv(
  mani,
  file.path(
    outdir,
    "01_analysis_manifest.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. Integrate TF results
# ============================================================

alltf <- do.call(
  rbind,
  integrated
)

rownames(
  alltf
) <- NULL


alltf$gene <- factor(
  alltf$gene,
  levels = rev(tfs)
)

alltf$dataset <- factor(
  alltf$dataset,
  levels = names(integrated)
)


write.csv(
  transform(
    alltf,
    gene = as.character(gene),
    dataset = as.character(dataset)
  ),
  file.path(
    outdir,
    "02_integrated_TF_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. Bubble heatmap
# ============================================================

plotdat <- alltf


plotdat$label <- ifelse(
  is.na(
    plotdat$P.Value
  ),
  "NA",
  plotdat$significance
)


# Clip extreme log2FC values for color visualization only
lim <- max(
  1,
  quantile(
    abs(
      plotdat$logFC
    ),
    0.95,
    na.rm = TRUE
  )
)


plotdat$plotFC <- pmax(
  -lim,
  pmin(
    lim,
    plotdat$logFC
  )
)


p <- ggplot(
  plotdat,
  aes(
    x = dataset,
    y = gene,
    fill = plotFC
  )
) +

  geom_point(
    shape = 21,
    size = 10,
    color = "#B8B8B8",
    stroke = 0.45
  ) +

  geom_text(
    aes(
      label = label
    ),
    size = 3.2,
    na.rm = TRUE
  ) +

  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    limits = c(
      -lim,
      lim
    ),
    name = "log2FC\n(color capped)"
  ) +

  labs(
    x = NULL,
    y = NULL,

    title = paste0(
      "Transcription-factor changes across ",
      "HCM, HF and TAC datasets"
    ),

    subtitle = paste0(
      "Direction: disease/TAC minus normal/sham; ",
      "* nominal P < 0.05, ** < 0.01, *** < 0.001"
    )
  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(
    panel.grid.major.x =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 50,
        hjust = 1,
        vjust = 1
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "right"
  )


# ============================================================
# 11. Save figure
# ============================================================

plot.width <- max(
  10,
  0.75 *
    length(
      unique(
        plotdat$dataset
      )
    ) +
    3
)


ggsave(
  file.path(
    outdir,
    "03_TF_bubble_heatmap.png"
  ),
  p,
  width = plot.width,
  height = 9,
  dpi = 320,
  bg = "white"
)


ggsave(
  file.path(
    outdir,
    "03_TF_bubble_heatmap.pdf"
  ),
  p,
  width = plot.width,
  height = 9,
  useDingbats = FALSE
)


# ============================================================
# 12. Analysis README
# ============================================================

writeLines(
  c(

    "MFAP5 transcription-factor cross-dataset analysis",

    paste(
      "Generated:",
      Sys.time()
    ),

    "",

    "Direction:",
    "Disease/TAC minus normal/sham.",

    "",

    "Methods:",
    paste0(
      "Count matrices: edgeR TMM normalization + ",
      "limma-voom."
    ),

    paste0(
      "Continuous matrices: limma after the ",
      "dataset-specific transformation indicated ",
      "in the analysis manifest."
    ),

    "",

    "Significance:",
    paste0(
      "Significance symbols in the bubble heatmap ",
      "use nominal P values."
    ),

    paste0(
      "Adjusted P values are retained in all ",
      "per-dataset and integrated result tables."
    ),

    "",

    "Dataset-specific notes:",

    paste0(
      "GSE133054 human HCM and HF were analyzed ",
      "as separate comparisons."
    ),

    paste0(
      "GSE133054 mouse TAC time points were pooled ",
      "and compared with sham samples."
    ),

    paste0(
      "GSE36961 and GSE249925 share many sample ",
      "identifiers and should not be treated as ",
      "independent cohorts in meta-analysis."
    ),

    "",

    "Visualization:",

    paste0(
      "Bubble colors represent log2 fold changes. ",
      "For visualization only, the color scale is ",
      "capped at the 95th percentile of absolute ",
      "log2FC values."
    ),

    paste0(
      "The exported CSV files retain the original ",
      "unmodified effect sizes."
    )
  ),

  file.path(
    outdir,
    "README_analysis.txt"
  )
)


# ============================================================
# 13. Completion message
# ============================================================

cat(
  "Completed ",
  nrow(mani),
  " comparisons and ",
  nrow(alltf),
  " TF rows.\n",
  sep = ""
)

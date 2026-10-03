# D --------------------
#!/usr/bin/env Rscript

# ============================================================
# Script: MFAP5_metabolite_classification_pie.R
#
# Description:
#   Classify differential metabolites from MFAP5 vs. Ctrl
#   into six selected KEGG metabolic categories and visualize
#   their composition using a pie chart.
#
# Metabolic categories:
#   1. Amino acid metabolism
#   2. Metabolism of other amino acids
#   3. Lipid metabolism
#   4. Nucleotide metabolism
#   5. Metabolism of cofactors and vitamins
#   6. Carbohydrate metabolism
#
# Input:
#   LC MFAP5-vs-Ctrl.diff.xls
#   3.xlsx
#
# Output:
#   1A20260914.pdf
#   1A_SuperClass_summary_20260914.csv
# ============================================================


# -----------------------------
# Load packages
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readxl)
})


# -----------------------------
# File paths
# -----------------------------
input.matches <- Sys.glob(
  "D:/Dataset/*/LC MFAP5-vs-Ctrl.diff.xls"
)

if (length(input.matches) != 1) {
  stop(
    "Expected exactly one LC MFAP5-vs-Ctrl.diff.xls under D:/Dataset"
  )
}

input.file <- input.matches[1]
base.dir <- dirname(input.file)

# "结果图"
figure.dir <- file.path(
  base.dir,
  intToUtf8(
    c(32467, 26524, 22270)
  )
)

mapping.file <- file.path(
  base.dir,
  "3.xlsx"
)

dir.create(
  figure.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# Read differential metabolites
# -----------------------------
dat <- read.delim(
  input.file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)


# -----------------------------
# Read KEGG classification
# -----------------------------
# Stage the Excel file at an ASCII-only temporary path
# to avoid potential path/encoding issues.

temp.mapping <- file.path(
  tempdir(),
  "metabolism_classification.xlsx"
)

if (
  !file.copy(
    mapping.file,
    temp.mapping,
    overwrite = TRUE
  )
) {
  stop(
    "Failed to stage 3.xlsx for reading."
  )
}

mapping <- read_xlsx(
  temp.mapping
)


# -----------------------------
# Check required columns
# -----------------------------
stopifnot(
  all(
    c(
      "Metabolites",
      "KEGG"
    ) %in% names(dat)
  ),
  all(
    c(
      "Classification_level2",
      "Substances"
    ) %in% names(mapping)
  )
)


# -----------------------------
# Selected metabolic classes
# -----------------------------
classes <- c(
  "Amino acid metabolism",
  "Metabolism of other amino acids",
  "Lipid metabolism",
  "Nucleotide metabolism",
  "Metabolism of cofactors and vitamins",
  "Carbohydrate metabolism"
)


# -----------------------------
# Build KEGG-to-class mapping
# -----------------------------
map.one <- mapping %>%
  filter(
    Classification_level2 %in% classes
  ) %>%
  select(
    Classification_level2,
    Substances
  ) %>%
  separate_rows(
    Substances,
    sep = ","
  ) %>%
  mutate(
    KEGG = str_trim(Substances),
    priority = match(
      Classification_level2,
      classes
    )
  ) %>%
  filter(
    !is.na(KEGG),
    KEGG != ""
  ) %>%
  arrange(priority) %>%
  distinct(
    KEGG,
    .keep_all = TRUE
  ) %>%
  select(
    KEGG,
    Classification_level2
  )


# -----------------------------
# Assign metabolites to classes
# -----------------------------
classified <- dat %>%
  filter(
    !is.na(KEGG),
    KEGG %in% map.one$KEGG
  ) %>%
  left_join(
    map.one,
    by = "KEGG"
  ) %>%
  distinct(
    Metabolites,
    .keep_all = TRUE
  )


# -----------------------------
# Summarize class composition
# -----------------------------
class.summary <- classified %>%
  count(
    Classification_level2,
    name = "Number"
  ) %>%
  complete(
    Classification_level2 = classes,
    fill = list(
      Number = 0
    )
  ) %>%
  mutate(
    Classification_level2 = factor(
      Classification_level2,
      levels = classes
    )
  ) %>%
  arrange(
    Classification_level2
  ) %>%
  mutate(
    Percent = Number / sum(Number) * 100,
    Label = paste0(
      as.character(Classification_level2),
      " (",
      Number,
      ", ",
      sprintf("%.1f", Percent),
      "%)"
    )
  )


# -----------------------------
# Colors
# -----------------------------
class.palette <- c(
  "#9467bd",
  "#17becf",
  "#aec7e8",
  "#ffbb78",
  "#ff9896",
  "#c49c94"
)

class.colors <- setNames(
  class.palette,
  class.summary$Label
)

class.summary$Label <- factor(
  class.summary$Label,
  levels = class.summary$Label
)


# -----------------------------
# Pie chart
# -----------------------------
p <- ggplot(
  class.summary,
  aes(
    x = "",
    y = Number,
    fill = Label
  )
) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.45
  ) +
  coord_polar(
    theta = "y",
    start = 0
  ) +
  scale_fill_manual(
    values = class.colors,
    drop = FALSE
  ) +
  labs(
    title = "Differential metabolites by metabolic class",
    subtitle = paste0(
      "MFAP5 vs. Ctrl (n = ",
      sum(class.summary$Number),
      ")"
    ),
    fill = "Metabolic class"
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      color = "grey35"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# -----------------------------
# Output paths
# -----------------------------
temp.pdf <- file.path(
  tempdir(),
  "1A20260914.pdf"
)

temp.csv <- file.path(
  tempdir(),
  "1A_SuperClass_summary_20260914.csv"
)


# -----------------------------
# Export figure
# -----------------------------
ggsave(
  temp.pdf,
  p,
  width = 11,
  height = 5.5,
  useDingbats = FALSE
)


# -----------------------------
# Export summary
# -----------------------------
write.csv(
  class.summary %>%
    mutate(
      Label = as.character(Label)
    ),
  temp.csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# -----------------------------
# Copy results to output folder
# -----------------------------
if (
  !file.copy(
    temp.pdf,
    file.path(
      figure.dir,
      "1A20260914.pdf"
    ),
    overwrite = TRUE
  )
) {
  stop(
    "Failed to copy PDF to output directory."
  )
}


if (
  !file.copy(
    temp.csv,
    file.path(
      figure.dir,
      "1A_SuperClass_summary_20260914.csv"
    ),
    overwrite = TRUE
  )
) {
  stop(
    "Failed to copy CSV to output directory."
  )
}


# -----------------------------
# Summary
# -----------------------------
message(
  "Completed: ",
  sum(class.summary$Number),
  paste0(
    " differential metabolites assigned uniquely ",
    "to 6 selected metabolic classes."
  )
)














# E --------------------
#!/usr/bin/env Rscript

# ============================================================
# Script: MFAP5_metabolomics_heatmap.R
# Description:
#   Identify significantly altered metabolites between MFAP5
#   and pcDNA groups and visualize:
#     1. KEGG metabolic classification
#     2. pcDNA metabolite expression
#     3. MFAP5 metabolite expression
#     4. Metabolite log2 fold changes
#
# Significance criteria:
#   P-value < 0.05
#   |log2FC| >= 0.585 (approximately FC >= 1.5 or <= 1/1.5)
#
# Input:
#   2.xlsx
#     Sheet1: differential metabolite statistics
#     Sheet2: metabolite abundance matrix
#   3.xlsx
#     Sheet1: KEGG metabolite classification
#
# Output:
#   PDF/PNG figure
#   Metabolites included in the figure
#   Plot-order table
# ============================================================


# -----------------------------
# Environment
# -----------------------------
invisible(
  try(
    Sys.setlocale("LC_CTYPE", "Chinese"),
    silent = TRUE
  )
)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})


# -----------------------------
# File paths
# -----------------------------
source.dir <- "D:/Dataset/代谢组学"

data.dir <- file.path(
  tempdir(),
  "metabolomics_plot_data"
)

output.dir <- file.path(
  data.dir,
  "results"
)

final.output.dir <- file.path(
  source.dir,
  "结果图"
)

dir.create(
  data.dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  output.dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  final.output.dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# -----------------------------
# Copy input files to ASCII path
# -----------------------------
ok <- file.copy(
  file.path(
    source.dir,
    c("2.xlsx", "3.xlsx")
  ),
  file.path(
    data.dir,
    c("2.xlsx", "3.xlsx")
  ),
  overwrite = TRUE
)

if (!all(ok)) {
  stop(
    "Failed to copy 2.xlsx or 3.xlsx to the temporary ASCII path."
  )
}


# -----------------------------
# Parameters
# -----------------------------
p.cutoff <- 0.05
logfc.cutoff <- 0.585

show.metabolite.names <- FALSE

date.tag <- "20260915"


# ============================================================
# 1. Differential metabolite analysis
# ============================================================

diff <- read_xlsx(
  file.path(data.dir, "2.xlsx"),
  sheet = "Sheet1"
) %>%
  distinct(
    Metabolites,
    .keep_all = TRUE
  ) %>%
  mutate(
    `-log10(q)` = -log10(
      pmax(
        `q-value`,
        .Machine$double.xmin
      )
    ),
    Significance = case_when(
      `p-value` < p.cutoff &
        log2FoldChange >= logfc.cutoff ~ "Significant Up",

      `p-value` < p.cutoff &
        log2FoldChange <= -logfc.cutoff ~ "Significant Down",

      TRUE ~ "Non-significant"
    )
  )


# ============================================================
# 2. Read metabolite abundance data
# ============================================================

expr.info <- read_xlsx(
  file.path(data.dir, "2.xlsx"),
  sheet = "Sheet2"
) %>%
  distinct(
    Metabolites,
    .keep_all = TRUE
  )


# ============================================================
# 3. Read KEGG metabolite classification
# ============================================================

path.raw <- read_xlsx(
  file.path(data.dir, "3.xlsx"),
  sheet = "Sheet1"
) %>%
  select(
    Classification_level2,
    Substances
  ) %>%
  separate_rows(
    Substances,
    sep = ","
  ) %>%
  mutate(
    across(
      c(
        Classification_level2,
        Substances
      ),
      str_trim
    )
  ) %>%
  filter(
    !is.na(Substances),
    Substances != "",
    !is.na(Classification_level2),
    Classification_level2 != ""
  ) %>%
  distinct()


# Preserve original KEGG class order
class.order <- unique(
  path.raw$Classification_level2
)


# Assign one metabolic class to each KEGG substance
path.one <- path.raw %>%
  mutate(
    class.rank = match(
      Classification_level2,
      class.order
    )
  ) %>%
  arrange(class.rank) %>%
  group_by(Substances) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    Substances,
    Classification_level2
  )


# ============================================================
# 4. Select significant metabolites with KEGG annotation
# ============================================================

plot.info <- expr.info %>%
  inner_join(
    diff %>%
      filter(
        Significance != "Non-significant"
      ) %>%
      select(
        Metabolites,
        log2FoldChange,
        `p-value`,
        `q-value`,
        `-log10(q)`,
        Significance
      ),
    by = "Metabolites"
  ) %>%
  inner_join(
    path.one,
    by = c(
      "KEGG" = "Substances"
    )
  ) %>%
  mutate(
    Classification_level2 = factor(
      Classification_level2,
      levels = class.order
    )
  ) %>%
  arrange(
    Classification_level2,
    desc(log2FoldChange),
    Metabolites
  ) %>%
  mutate(
    RowID = make.unique(
      as.character(Metabolites)
    ),
    Metabolite_plot = factor(
      RowID,
      levels = rev(RowID)
    )
  )


if (nrow(plot.info) == 0) {
  stop(
    paste0(
      "No significant metabolites matched both ",
      "the expression table and KEGG classification."
    )
  )
}


# ============================================================
# 5. Identify sample columns
# ============================================================

mfap.cols <- grep(
  "^MFAP5_",
  names(plot.info),
  value = TRUE
)

ctrl.cols <- grep(
  "^pcdna",
  names(plot.info),
  value = TRUE,
  ignore.case = TRUE
)

if (
  length(mfap.cols) == 0 ||
  length(ctrl.cols) == 0
) {
  stop(
    paste0(
      "MFAP5 or pcDNA sample columns were ",
      "not detected in 2.xlsx Sheet2."
    )
  )
}


# ============================================================
# 6. Row-wise Z-score normalization
# ============================================================

expr.mat <- as.matrix(
  plot.info[
    ,
    c(
      ctrl.cols,
      mfap.cols
    )
  ]
)

storage.mode(expr.mat) <- "numeric"


expr.scale <- t(
  apply(
    expr.mat,
    1,
    function(x) {

      s <- sd(
        x,
        na.rm = TRUE
      )

      if (
        !is.finite(s) ||
        s == 0
      ) {
        rep(
          0,
          length(x)
        )
      } else {
        (
          x -
            mean(
              x,
              na.rm = TRUE
            )
        ) / s
      }
    }
  )
)


colnames(expr.scale) <- colnames(expr.mat)
rownames(expr.scale) <- plot.info$RowID


# Convert expression matrix to long format
heat.long <- as.data.frame(
  expr.scale
) %>%
  mutate(
    RowID = rownames(expr.scale)
  ) %>%
  pivot_longer(
    -RowID,
    names_to = "Sample",
    values_to = "Zscore"
  ) %>%
  left_join(
    plot.info %>%
      select(
        RowID,
        Metabolite_plot,
        Classification_level2
      ),
    by = "RowID"
  ) %>%
  mutate(
    Group = if_else(
      Sample %in% mfap.cols,
      "MFAP5",
      "pcDNA"
    ),
    Sample = factor(
      Sample,
      levels = c(
        ctrl.cols,
        mfap.cols
      )
    )
  )


# ============================================================
# 7. KEGG class settings
# ============================================================

plot.info$Classification_level2 <- droplevels(
  plot.info$Classification_level2
)

heat.long$Classification_level2 <- factor(
  heat.long$Classification_level2,
  levels = levels(
    plot.info$Classification_level2
  )
)

classes <- levels(
  plot.info$Classification_level2
)

class.palette <- c(
  "#9467bd",
  "#17becf",
  "#aec7e8",
  "#ffbb78",
  "#ff9896",
  "#c49c94"
)

class.colors <- setNames(
  rep(
    class.palette,
    length.out = length(classes)
  ),
  classes
)


# Select middle metabolite for class label
label.data <- plot.info %>%
  group_by(
    Classification_level2
  ) %>%
  slice(
    ceiling(n() / 2)
  ) %>%
  ungroup()


# ============================================================
# 8. Shared plot settings
# ============================================================

facet.spec <- facet_grid(
  rows = vars(
    Classification_level2
  ),
  scales = "free_y",
  space = "free_y"
)


theme.panel <- theme_minimal(
  base_size = 9
) +
  theme(
    panel.grid = element_blank(),

    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.55
    ),

    panel.spacing.y = unit(
      1.5,
      "mm"
    ),

    strip.text.y = element_blank(),
    strip.background = element_blank(),

    axis.title.y = element_blank(),

    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 11
    ),

    plot.margin = margin(
      3,
      3,
      3,
      3
    )
  )


# ============================================================
# 9. Metabolic classification panel
# ============================================================

p.class <- ggplot(
  plot.info,
  aes(
    x = 1,
    y = Metabolite_plot,
    fill = Classification_level2
  )
) +
  geom_tile(
    width = 1,
    height = 1
  ) +
  geom_text(
    data = label.data,
    aes(
      label = Classification_level2
    ),
    size = 2.6,
    lineheight = 0.9
  ) +
  facet.spec +
  scale_fill_manual(
    values = class.colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  labs(
    title = "Metabolism class",
    x = NULL
  ) +
  theme.panel +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),

    axis.text.y = element_text(
      size = 5.4,
      color = "black",
      margin = margin(r = 3)
    ),

    legend.position = "none"
  )


# ============================================================
# 10. Expression heatmap function
# ============================================================

make.heat <- function(
  group.name,
  title.name
) {

  ggplot(
    filter(
      heat.long,
      Group == group.name
    ),
    aes(
      x = Sample,
      y = Metabolite_plot,
      fill = Zscore
    )
  ) +
    geom_tile() +
    facet.spec +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -max(
          abs(heat.long$Zscore),
          na.rm = TRUE
        ),
        max(
          abs(heat.long$Zscore),
          na.rm = TRUE
        )
      ),
      oob = squish,
      name = "Z-score"
    ) +
    labs(
      title = title.name,
      x = NULL
    ) +
    theme.panel +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),

      axis.text.y = if (
        show.metabolite.names
      ) {
        element_text(size = 5)
      } else {
        element_blank()
      },

      axis.ticks.y = element_blank(),
      legend.position = "top"
    )
}


p.ctrl <- make.heat(
  "pcDNA",
  "pcDNA"
)

p.mfap <- make.heat(
  "MFAP5",
  "MFAP5"
)


# ============================================================
# 11. log2FC panel
# ============================================================

x.limit <- max(
  abs(plot.info$log2FoldChange),
  na.rm = TRUE
) * 1.08


p.fc <- ggplot(
  plot.info,
  aes(
    y = Metabolite_plot
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey65",
    linewidth = 0.45
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      yend = Metabolite_plot
    ),
    linetype = "dotted",
    colour = "grey45",
    linewidth = 0.45
  ) +
  geom_point(
    aes(
      x = log2FoldChange,
      colour = `-log10(q)`,
      size = abs(log2FoldChange)
    ),
    alpha = 0.95
  ) +
  facet.spec +
  scale_x_continuous(
    limits = c(
      -x.limit,
      x.limit
    ),
    breaks = pretty(
      c(
        -x.limit,
        x.limit
      ),
      n = 5
    )
  ) +
  scale_colour_gradient(
    low = "#F7B6D2",
    high = "#C40064",
    name = "-log10(q)"
  ) +
  scale_size_continuous(
    range = c(
      1.5,
      4.2
    ),
    name = "|Log2FC|"
  ) +
  labs(
    title = "Log2FC",
    x = "Log2FC"
  ) +
  theme.panel +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )


# ============================================================
# 12. Combine panels
# ============================================================

final.plot <-
  p.class +
  p.ctrl +
  p.mfap +
  p.fc +
  plot_layout(
    widths = c(
      3.8,
      2.1,
      2.1,
      2.3
    ),
    guides = "collect"
  )


# Dynamic figure height
fig.height <- max(
  7,
  2.5 +
    nrow(plot.info) * 0.115
)


# ============================================================
# 13. Export figures
# ============================================================

ggsave(
  file.path(
    output.dir,
    paste0(
      "1E2_metabolism_class_heatmap_Log2FC_",
      "metabolite_names_",
      date.tag,
      ".pdf"
    )
  ),
  final.plot,
  width = 15.5,
  height = fig.height,
  limitsize = FALSE
)


ggsave(
  file.path(
    output.dir,
    paste0(
      "1E2_metabolism_class_heatmap_Log2FC_",
      "metabolite_names_",
      date.tag,
      ".png"
    )
  ),
  final.plot,
  width = 15.5,
  height = fig.height,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)


# ============================================================
# 14. Export metabolite information
# ============================================================

write.csv(
  plot.info %>%
    select(
      Metabolites,
      HMDB,
      KEGG,
      Classification_level2,
      log2FoldChange,
      `p-value`,
      `q-value`,
      Significance
    ),
  file.path(
    output.dir,
    paste0(
      "1E2_metabolites_used_",
      date.tag,
      ".csv"
    )
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 15. Export metabolite plot order
# ============================================================

plot.order.txt <- plot.info %>%
  mutate(
    Plot_Order = row_number()
  ) %>%
  select(
    Plot_Order,
    Metabolites,
    HMDB,
    KEGG,
    Classification_level2,
    log2FoldChange,
    `p-value`,
    `q-value`,
    Significance
  )


temp.txt <- file.path(
  data.dir,
  paste0(
    "1E2_metabolites_plot_order_",
    date.tag,
    ".txt"
  )
)


write.table(
  plot.order.txt,
  temp.txt,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  fileEncoding = "UTF-8",
  na = ""
)


if (
  !file.copy(
    temp.txt,
    file.path(
      source.dir,
      basename(temp.txt)
    ),
    overwrite = TRUE
  )
) {
  warning(
    "The plot-order TXT file could not be copied to the source directory."
  )
}


# ============================================================
# 16. Copy results to final output directory
# ============================================================

result.files <- list.files(
  output.dir,
  full.names = TRUE
)

if (
  !all(
    file.copy(
      result.files,
      file.path(
        final.output.dir,
        basename(result.files)
      ),
      overwrite = TRUE
    )
  )
) {
  warning(
    "Some result files could not be copied to the final directory."
  )
}


# ============================================================
# Summary
# ============================================================

message(
  "Completed: ",
  nrow(plot.info),
  " significant metabolites with KEGG classifications were included."
)













# G ----------------------------
#!/usr/bin/env Rscript

# ============================================================
# Script: MFAP5_energy_metabolites_violin_boxplot.R
#
# Description:
#   Compare selected energy-related metabolites between
#   MFAP5 and control groups using violin/box plots.
#
# Metabolites:
#   - Lactate
#   - ATP
#   - Pyruvic acid
#
# Statistical test:
#   Two-sided Welch's t-test
#
# Input:
#   MFAP5_vs_ctrl_heatmap_row_cluster.xlsx
#
# Output:
#   energy_metabolites_violin_boxplot_20260915.pdf
#   energy_metabolites_violin_boxplot_20260915.png
#   energy_metabolites_violin_boxplot_stats_20260915.txt
# ============================================================


# -----------------------------
# Environment
# -----------------------------
invisible(
  try(
    Sys.setlocale("LC_CTYPE", "Chinese"),
    silent = TRUE
  )
)

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})


# -----------------------------
# Locate input file
# -----------------------------
input.matches <- Sys.glob(
  paste0(
    "D:/Dataset/*/MWY-25-8715-a12*/MWY-25-8715-a/",
    "2.Basic_Analysis/Difference_analysis/",
    "MFAP5_vs_ctrl/heatmap/",
    "MFAP5_vs_ctrl_heatmap_row_cluster.xlsx"
  )
)

if (length(input.matches) != 1) {
  stop(
    "Expected exactly one MFAP5_vs_ctrl heatmap workbook."
  )
}

input.file <- input.matches[1]


# -----------------------------
# Define project/output directory
# -----------------------------
metabolism.dir <- input.file

for (i in seq_len(7)) {
  metabolism.dir <- dirname(metabolism.dir)
}

# "结果图"
figure.dir <- file.path(
  metabolism.dir,
  intToUtf8(
    c(32467, 26524, 22270)
  )
)

dir.create(
  figure.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# Stage Excel file to ASCII path
# -----------------------------
temp.xlsx <- file.path(
  tempdir(),
  "MFAP5_vs_ctrl_heatmap_row_cluster.xlsx"
)

if (
  !file.copy(
    input.file,
    temp.xlsx,
    overwrite = TRUE
  )
) {
  stop(
    "Failed to stage input workbook."
  )
}


# -----------------------------
# Read data
# -----------------------------
dat <- read_xlsx(
  temp.xlsx,
  sheet = "sheet1"
)


# -----------------------------
# Sample and metabolite settings
# -----------------------------
sample.cols <- c(
  paste0("ctrl-", 1:6),
  paste0("MFAP5-", 1:6)
)

targets <- c(
  "Lactate",
  "ATP",
  "Pyruvic acid"
)


# -----------------------------
# Check required data
# -----------------------------
stopifnot(
  all(
    c(
      "Compounds",
      sample.cols
    ) %in% names(dat)
  ),
  all(
    targets %in% dat$Compounds
  )
)


# -----------------------------
# Prepare long-format data
# -----------------------------
plot.data <- dat %>%
  filter(
    Compounds %in% targets
  ) %>%
  select(
    Compounds,
    all_of(sample.cols)
  ) %>%
  pivot_longer(
    -Compounds,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  mutate(
    Abundance = as.numeric(Abundance),

    Group = if_else(
      grepl(
        "^MFAP5-",
        Sample
      ),
      "MFAP5",
      "ctrl"
    ),

    Group = factor(
      Group,
      levels = c(
        "ctrl",
        "MFAP5"
      )
    ),

    Compounds = factor(
      Compounds,
      levels = targets
    )
  )


if (anyNA(plot.data$Abundance)) {
  stop(
    "Missing or nonnumeric abundance values detected."
  )
}


# ============================================================
# Statistical analysis
# ============================================================

stats <- plot.data %>%
  group_by(
    Compounds
  ) %>%
  summarise(
    n.MFAP5 = sum(
      Group == "MFAP5"
    ),

    n.ctrl = sum(
      Group == "ctrl"
    ),

    mean.MFAP5 = mean(
      Abundance[
        Group == "MFAP5"
      ]
    ),

    mean.ctrl = mean(
      Abundance[
        Group == "ctrl"
      ]
    ),

    difference =
      mean.MFAP5 -
      mean.ctrl,

    p.value = t.test(
      Abundance[
        Group == "MFAP5"
      ],
      Abundance[
        Group == "ctrl"
      ],
      var.equal = FALSE
    )$p.value,

    .groups = "drop"
  ) %>%
  mutate(
    stars = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),

    p.text = if_else(
      p.value < 0.001,
      formatC(
        p.value,
        format = "e",
        digits = 2
      ),
      formatC(
        p.value,
        format = "f",
        digits = 4
      )
    ),

    label = paste0(
      "p = ",
      p.text,
      " (",
      stars,
      ")"
    )
  )


# -----------------------------
# Determine annotation positions
# -----------------------------
ranges <- plot.data %>%
  group_by(
    Compounds
  ) %>%
  summarise(
    y.min = min(
      Abundance
    ),

    y.max = max(
      Abundance
    ),

    .groups = "drop"
  ) %>%
  mutate(
    span = pmax(
      y.max - y.min,
      0.25
    ),

    y.bracket =
      y.max +
      0.12 * span,

    y.text =
      y.max +
      0.22 * span
  )


stats <- left_join(
  stats,
  ranges,
  by = "Compounds"
)


# ============================================================
# Plot settings
# ============================================================

group.colors <- c(
  MFAP5 = "#1B9E77",
  ctrl = "#D95F02"
)


# ============================================================
# Violin + boxplot + individual samples
# ============================================================

p <- ggplot(
  plot.data,
  aes(
    x = Group,
    y = Abundance,
    fill = Group,
    color = Group
  )
) +

  # Violin distribution
  geom_violin(
    width = 0.68,
    trim = FALSE,
    alpha = 0.48,
    linewidth = 0.55
  ) +

  # Boxplot
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.78,
    color = "black",
    linewidth = 0.5
  ) +

  # Individual samples
  geom_point(
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 20260915
    ),
    size = 1.5,
    alpha = 0.9
  ) +

  # Significance bracket
  geom_segment(
    data = stats,
    aes(
      x = 1,
      xend = 2,
      y = y.bracket,
      yend = y.bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +

  geom_segment(
    data = stats,
    aes(
      x = 1,
      xend = 1,
      y = y.bracket -
        0.04 * span,
      yend = y.bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +

  geom_segment(
    data = stats,
    aes(
      x = 2,
      xend = 2,
      y = y.bracket -
        0.04 * span,
      yend = y.bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +

  # P-value
  geom_text(
    data = stats,
    aes(
      x = 1.5,
      y = y.text,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +

  # One panel per metabolite
  facet_wrap(
    ~Compounds,
    nrow = 1,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = group.colors
  ) +

  scale_color_manual(
    values = group.colors
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.06,
        0.22
      )
    )
  ) +

  labs(
    x = NULL,
    y = "Standardized abundance",
    title = "Selected energy metabolites"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    legend.position = "none",

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face = "bold",
        size = 11
      ),

    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),

    axis.text.x =
      element_text(
        face = "bold"
      ),

    panel.spacing =
      grid::unit(
        0.85,
        "lines"
      )
  )


# ============================================================
# Export figure
# ============================================================

temp.pdf <- file.path(
  tempdir(),
  "energy_metabolites_violin_boxplot_20260915.pdf"
)

temp.png <- file.path(
  tempdir(),
  "energy_metabolites_violin_boxplot_20260915.png"
)


ggsave(
  temp.pdf,
  p,
  width = 7.6,
  height = 3.8,
  useDingbats = FALSE
)

ggsave(
  temp.png,
  p,
  width = 7.6,
  height = 3.8,
  dpi = 400
)


# -----------------------------
# Final output files
# -----------------------------
final.pdf <- file.path(
  figure.dir,
  "energy_metabolites_violin_boxplot_20260915.pdf"
)

final.png <- file.path(
  figure.dir,
  "energy_metabolites_violin_boxplot_20260915.png"
)

final.stats <- file.path(
  figure.dir,
  "energy_metabolites_violin_boxplot_stats_20260915.txt"
)


if (
  !file.copy(
    temp.pdf,
    final.pdf,
    overwrite = TRUE
  ) ||
  !file.copy(
    temp.png,
    final.png,
    overwrite = TRUE
  )
) {
  stop(
    "Failed to copy plot outputs."
  )
}


# ============================================================
# Export statistics
# ============================================================

write.table(
  stats %>%
    select(
      Compounds,
      n.MFAP5,
      n.ctrl,
      mean.MFAP5,
      mean.ctrl,
      difference,
      p.value,
      stars
    ),
  final.stats,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# -----------------------------
# Print summary
# -----------------------------
print(
  stats %>%
    select(
      Compounds,
      n.MFAP5,
      n.ctrl,
      mean.MFAP5,
      mean.ctrl,
      difference,
      p.value,
      stars
    )
)

message(
  "Saved plot and statistics to: ",
  figure.dir
)

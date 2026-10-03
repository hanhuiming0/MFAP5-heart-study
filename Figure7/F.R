# F ----------------
#!/usr/bin/env Rscript

# ============================================================
# Script: SCP1303_KLF12_positive_celltype_composition.R
#
# Description:
#   Identify KLF12-positive cells in the SCP1303 single-cell/
#   single-nucleus RNA-seq dataset and characterize their
#   cell-type composition in Normal and HCM samples.
#
# Definition:
#   KLF12-positive cells: KLF12 expression > 0
#
# Analyses:
#   1. Extract KLF12-positive cells
#   2. Calculate cell-type composition within KLF12+ cells
#      separately for Normal and HCM
#   3. Visualize KLF12+ cells on UMAP
#   4. Visualize cell-type proportions using stacked bar plots
#
# Output:
#   - KLF12_positive.Sample_Type.CellType.Normal_HCM.pdf
#   - KLF12_positive_CellType_proportion_Normal_HCM.pdf
#   - KLF12_positive_CellType_proportion_Normal_HCM.csv
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
})


# ============================================================
# 2. File paths
# ============================================================

input.file <- paste0(
  "/storage/data/Hhm/Heart_Study/data/sc_sn/",
  "SCP1303/Result/SCP1303.Normal_HCM.Cell_Type.RDS"
)

output.dir <- "."

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Parameters
# ============================================================

target.gene <- "KLF12"

sample.types <- c(
  "Normal",
  "HCM"
)


# ============================================================
# 4. Cell-type colors
# ============================================================

CellType_Color <- c(
  "Adipocyte"     = "#E69F00",
  "Cardiomyocyte" = "#D55E00",
  "Endothelial"   = "#56B4E9",
  "Fibroblast"    = "#009E73",
  "Mast"          = "#CC79A7",
  "Myeloid"       = "#0072B2",
  "Neuronal"      = "#F0E442",
  "Pericyte"      = "#9467BD",
  "Smooth muscle" = "#8C564B",
  "T/NK"          = "#17BECF"
)


# ============================================================
# 5. Load Seurat object
# ============================================================

SCP1303 <- readRDS(
  input.file
)

DefaultAssay(SCP1303) <- "RNA"


# Check required metadata
required.metadata <- c(
  "Sample_Type",
  "Cell_Type"
)

if (
  !all(
    required.metadata %in%
      colnames(SCP1303@meta.data)
  )
) {
  stop(
    "Sample_Type or Cell_Type was not found in the metadata."
  )
}


# Check target gene
if (
  !target.gene %in%
    rownames(SCP1303)
) {
  stop(
    paste0(
      target.gene,
      " was not found in the RNA assay."
    )
  )
}


# ============================================================
# 6. Extract KLF12-positive cells
# ============================================================

SCP1303.KLF12.pos <- subset(
  SCP1303,
  subset = KLF12 > 0
)


# Define sample order
SCP1303.KLF12.pos$Sample_Type <- factor(
  SCP1303.KLF12.pos$Sample_Type,
  levels = sample.types
)


message(
  "Total KLF12+ cells: ",
  ncol(SCP1303.KLF12.pos)
)


# ============================================================
# 7. Extract metadata of KLF12-positive cells
# ============================================================

KLF12.pos.meta <- SCP1303.KLF12.pos@meta.data

KLF12.pos.meta$Sample_Type <- factor(
  KLF12.pos.meta$Sample_Type,
  levels = sample.types
)


# ============================================================
# 8. Calculate cell-type composition
# ============================================================

KLF12.pos.ratio <- KLF12.pos.meta %>%
  filter(
    Sample_Type %in% sample.types,
    !is.na(Cell_Type)
  ) %>%
  count(
    Sample_Type,
    Cell_Type,
    name = "Cell_Number"
  ) %>%
  group_by(
    Sample_Type
  ) %>%
  mutate(
    Total_KLF12_Positive = sum(
      Cell_Number
    ),
    Proportion =
      Cell_Number /
      Total_KLF12_Positive
  ) %>%
  ungroup()


print(
  KLF12.pos.ratio
)


# ============================================================
# 9. UMAP visualization of KLF12-positive cells
# ============================================================

p.umap <- DimPlot(
  SCP1303.KLF12.pos,
  group.by = "Cell_Type",
  cols = CellType_Color,
  split.by = "Sample_Type"
)


ggsave(
  filename = file.path(
    output.dir,
    "KLF12_positive.Sample_Type.CellType.Normal_HCM.pdf"
  ),
  plot = p.umap,
  width = 6,
  height = 3.5,
  useDingbats = FALSE
)


# ============================================================
# 10. Cell-type composition plot
# ============================================================

p1 <- ggplot(
  KLF12.pos.ratio,
  aes(
    x = Sample_Type,
    y = Proportion,
    fill = Cell_Type
  )
) +

  geom_col(
    width = 0.7,
    color = "white",
    linewidth = 0.3
  ) +

  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  scale_fill_manual(
    values = CellType_Color,
    na.value = "grey70"
  ) +

  labs(
    title = "Cell-type composition of KLF12+ cells",
    x = NULL,
    y = "Proportion of KLF12+ cells",
    fill = "Cell type"
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),

    axis.text.x = element_text(
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    legend.title = element_text(
      face = "bold"
    )
  )


print(p1)


# ============================================================
# 11. Save cell-type composition plot
# ============================================================

ggsave(
  filename = file.path(
    output.dir,
    "KLF12_positive_CellType_proportion_Normal_HCM.pdf"
  ),
  plot = p1,
  width = 4,
  height = 4,
  useDingbats = FALSE
)


# ============================================================
# 12. Export cell-type composition table
# ============================================================

write.csv(
  KLF12.pos.ratio,
  file = file.path(
    output.dir,
    "KLF12_positive_CellType_proportion_Normal_HCM.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. Summary
# ============================================================

message(
  "Cell types detected in KLF12+ cells: ",
  paste(
    sort(
      unique(
        SCP1303.KLF12.pos$Cell_Type
      )
    ),
    collapse = ", "
  )
)

message(
  "Analysis completed."
)

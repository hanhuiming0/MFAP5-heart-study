library(CellChat)
## C -----------
figshare_5777948 <- readRDS('/storage/data/Hhm/Heart_Study/data/sc_sn/5.figshare_5777948/Result/figshare_5777948_Cell_Type.RDS')
# ①.按照Sample_Type分组分析 
sc.data <- figshare_5777948 
{
  sc.data$Cell_Type <- as.character(sc.data$Cell_Type)
  
  rna.data <- GetAssayData(
    sc.data,
    assay = "RNA",
    layer = "data"
  )
  
  fib.cells <- colnames(sc.data)[
    sc.data$Cell_Type == "Fibroblast"
  ]
  
  MFAP5.pos.cells <- fib.cells[
    rna.data["MFAP5", fib.cells] > 0
  ]
  
  sc.data$Cell_Type[
    sc.data$Cell_Type == "Fibroblast"
  ] <- "MFAP5- Fibroblast"
  
  sc.data$Cell_Type[
    colnames(sc.data) %in% MFAP5.pos.cells
  ] <- "MFAP5+ Fibroblast"
  
  table(sc.data$Sample_Type, sc.data$Cell_Type)
  
  # ②. 按 Sample_Type 拆分
  sample.groups <- unique(as.character(sc.data$Sample_Type))
  
  cellchat.list <- setNames(
    vector("list", length(sample.groups)),
    sample.groups
  )
  
  for (group.name in sample.groups) {
    message("Running CellChat: ", group.name)
    cells.use <- colnames(sc.data)[
      sc.data$Sample_Type == group.name
    ]
    object.use <- subset(
      sc.data,
      cells = cells.use
    )
    # 表达矩阵
    expr.use <- GetAssayData(
      object.use,
      assay = "RNA",
      layer = "data"
    )
    # CellChat metadata
    meta.use <- data.frame(
      labels = as.character(object.use$Cell_Type),
      row.names = colnames(object.use)
    )
    cellchat.use <- createCellChat(
      object = expr.use,
      meta = meta.use,
      group.by = "labels"
    )
    cellchat.use@DB <- CellChatDB.human
    
    cellchat.use <- subsetData(cellchat.use)
    cellchat.use <- identifyOverExpressedGenes(cellchat.use)
    cellchat.use <- identifyOverExpressedInteractions(cellchat.use)
    
    cellchat.use <- computeCommunProb(
      cellchat.use,
      type = "triMean"
    )
    
    cellchat.use <- filterCommunication(
      cellchat.use,
      min.cells = 10
    )
    
    cellchat.use <- computeCommunProbPathway(cellchat.use)
    cellchat.use <- aggregateNet(cellchat.use)
    
    cellchat.list[[group.name]] <- cellchat.use
  }
  
  # ③. 合并两组 CellChat 对象
  cellchat.merged <- mergeCellChat(
    cellchat.list,
    add.names = names(cellchat.list),
    cell.prefix = TRUE
  )
  saveRDS(
    cellchat.list,
    file = "/storage/data/Hhm/Heart_Study/data/sc_sn/5.figshare_5777948/results/figshare_5777948.cellchat_by_SampleType.RDS"
  )
  saveRDS(
    cellchat.merged,
    file = "/storage/data/Hhm/Heart_Study/data/sc_sn/5.figshare_5777948/results/figshare_5777948.cellchat_merged.RDS"
  )

comparison.use <- match(
  c("Normal", "HCM"),
  names(cellchat.merged@net)
)
focus.cell <- "MFAP5+ Fibroblast"
cell.types <- levels(cellchat.list[["Normal"]]@idents)
other.cells <- setdiff(cell.types, focus.cell)
pdf("/storage/data/Hhm/Heart_Study/data/sc_sn/5.figshare_5777948/figures/HCM_vs_Normal_diffInteraction_count.pdf",
  width = 5.3,
  height = 6.4
)
netVisual_diffInteraction(
  cellchat.merged,
  comparison = comparison.use,
  measure = "count",
  sources.use = focus.cell,
  targets.use = other.cells,
  weight.scale = TRUE,
  label.edge = TRUE,
  title.name = "MFAP5+ Fibroblast outgoing: HCM vs. Normal")
dev.off()

# D -------------------------
pdf("/storage/data/Hhm/Heart_Study/data/sc_sn/5.figshare_5777948/figures/figshare_5777948.CellChat_netVisual_heatmap.pdf",width = 4.5,height = 4.2)
netVisual_heatmap(
  cellchat.merged,
  comparison = c(1, 2),
  measure = "count",
  color.use = CellType_Color,
  color.heatmap = c("#f0f7fc", "#4a88b7")
)
netVisual_heatmap(
  cellchat.merged,
  comparison = c(1, 2),
  measure = "weight",
  color.use = CellType_Color,
  color.heatmap = c("#f0f7fc", "#58a996")
)
dev.off()

library(Seurat)
setwd('/storage/data/Hhm/Heart_Study/data/sc_sn/SCP1303')
SCP1303 <- readRDS('/storage/data/Hhm/Heart_Study/data/sc_sn/SCP1303/Result/SCP1303_Cell_Type.RDS')
## A-----------
CellType_Color <- c('Cardiomyocyte'="#D42727",'Fibroblast'="#2078B3",
                    'Endothelial'="#279E68",'Myeloid'="#FCBD76",
                    'Pericyte = '#458A74', 'Smooth muscle' = '#91C6C2',
                    'T/NK'="#12C0D0", 'Neuronal' = "#AFC6E7",
                    'Mast'="#8D554F",'Adipocyte'="#FF7D0E")
pdf('/storage/data/Hhm/Heart_Study/NZH/SCP1303.celltype.umap.pdf', width = 3.8,height = 3.8)
DimPlot(SCP1303, reduction = "umap",group.by = 'CellType',
        cols = CellType_Color, label = T) + NoLegend()
dev.off()
pdf('/storage/data/Hhm/Heart_Study/NZH/SCP1303.MFAP5.umap.pdf', width = 3.8,height = 3.8)
FeaturePlot(SCP1303, reduction = "umap",features = 'MFAP5',
            cols = c("grey", "#7301A8FF")) + NoLegend()
dev.off()

## B -----------
pdf('/storage/data/Hhm/Heart_Study/NZH/SCP1303.MFAP5.DotPlot.umap.pdf',
    width = 4,height = 3.8)
DotPlot(SCP1303,features = 'MFAP5',cols = c("grey", "#7301A8FF"),
        group.by = 'Cell_Type')
dev.off()


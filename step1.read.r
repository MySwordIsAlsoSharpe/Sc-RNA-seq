## 1. read h5 ##
h5     = paste0('rawdata/', names(samples), '/outs/filtered_feature_bc_matrix.h5')
counts = do.call(cbind, lapply(seq(h5), function(i) {
    f    = h5[i]
    name = samples[i]
    message('Read ', name, ': ', f)
    count = Read10X_h5(f) 
    colnames(count) = paste0(name, '_', colnames(count))
    count
}))
obj    = CreateSeuratObject(counts); rm(counts)
#### cell cycle
cc.s   = intersect(c(str_to_title(cc.genes.updated.2019$s.genes),   cc.genes.updated.2019$s.genes),   rownames(obj))
cc.g2m = intersect(c(str_to_title(cc.genes.updated.2019$g2m.genes), cc.genes.updated.2019$g2m.genes), rownames(obj))
obj    = NormalizeData(obj)
obj    = CellCycleScoring(obj, s.features = cc.s, g2m.features = cc.g2m)
obj$cc.diff = obj$S.Score - obj$G2M.Score
#### subset
obj = subset(obj, nFeature_RNA > 200 & nCount_RNA > 1000 & nCount_RNA < 70000 & mt.pct < 20)

## 2. Normalize ##
obj = Seurat.Normalize(obj)
obj = Seurat.Cluster(obj, group.by = 'samples')
p   = DimPlot(obj, label = T) + theme(text = element_text(family = 'serif')); p

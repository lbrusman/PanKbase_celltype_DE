suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))

outdir <- "/tscc/nfs/home/lebrusman/Gaulton_lab/code/pankbase_celltype_marker_genes/matrices/"

data <- readRDS("/tscc/nfs/home/lebrusman/Gaulton_lab/data/ADA_object/250424_ADA_object_metadata_v3_3.rds")
print(data)

for (s in c("SRR27326986", "SRR27326987", "SRR27326992", "SRR27326993",
            "SRR27326994", "SRR27326995", "SRR27326996", "SRR27326997")) {
    data@meta.data[data@meta.data$samples == s, "samples"] <- paste0(data@meta.data[data@meta.data$samples == s, "samples"],
                                                                    "__", data@meta.data[data@meta.data$samples == s, "treatments"])
}

data@meta.data$coarse_annot <- gsub(" ", "", data@meta.data$coarse_annot)
Idents(data) <- data@meta.data$coarse_annot
cell.types <- unique(data$coarse_annot)

Idents(data) <- data@meta.data$coarse_annot

for (c in cell.types) {
    print(c)
    c_seurat <- subset(data, idents = c)
    allbut_seurat <- subset(data, idents = c, invert = TRUE)

    c_bulk <- AggregateExpression(c_seurat, group.by = "samples", return.seurat = TRUE)[["RNA"]]$counts %>% as.data.frame()
    c_bulk <- cbind(gene = rownames(c_bulk), c_bulk)
    allbut_bulk <- AggregateExpression(allbut_seurat, group.by = "samples", return.seurat = TRUE)[["RNA"]]$counts %>% as.data.frame()
    allbut_bulk <- cbind(gene = rownames(allbut_bulk), allbut_bulk)

    write.table(c_bulk, paste0(outdir, "cell_mtx/", c, "_persample_RNA_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(allbut_bulk, paste0(outdir, "allbut_mtx/", c, "_persample_RNA_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
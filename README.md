# PanKbase_celltype_DE
This repo houses the between-cell type differential expression analysis for PanKbase.

This code uses DESeq2 to perform differential expression analysis between sample-level pseudobulk matrices for the cell type of interest vs. all other cells. It then uses a ranked list of all genes (after removal of mitochondrial and ribosomal genes) to run fGSEA (https://bioconductor.org/packages/release/bioc/html/fgsea.html).

## General filtering

Donors were filtered based on the following criteria:
- Must be non-diabetic and autoantibody negative (if there is autoantibody data available)
- Samples must not be treated with any additional reagents (`treatments == "no_treatment"`)
- Must have $\geq$ 20 cells in the cell type of interest

Genes were filtered based on the following criteria:
- Must have $\geq$ 5 counts in $\geq$ 25% of samples

## Contrast info
DESeq formula:
- `~ subtype + sample` (where `subtype` is cell type of `"interest"` vs. `"other"`)

## Output files
- `<celltype>_marker_genes_deseq.tsv`: DESeq results comparing cell type of interest to all other cell types. Positive log2FoldChange indicates a gene is more highly expressed in that cell type compared to other cells.
- `<celltype>_fGSEA_res_all.tsv`: Results from fGSEA for *all* pathways.
- `<celltype>_fGSEA_res_signif.tsv`: Results from fGSEA for *significant* pathways (p < 0.1).
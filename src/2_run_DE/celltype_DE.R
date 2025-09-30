suppressMessages(library(plyr))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))
suppressMessages(library(knitr))
suppressMessages(library(stringr))
suppressMessages(library(DESeq2))
suppressPackageStartupMessages(library(fgsea))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
# opts_chunk$set(tidy=TRUE)
warnLevel <- getOption('warn')
options(warn = -1)

#input dir
dir <- "/tscc/nfs/home/lebrusman/Gaulton_lab/code/pankbase_celltype_marker_genes/matrices/"
files <- list.files(dir, pattern="_persample_RNA_counts.tsv")
cells <- gsub("_persample_RNA_counts.tsv", "", files)

#get cell name columns
#for pankbase
cell_names <- gsub("\\+", ".", cells)
cell_names <- gsub("\\(", ".", cell_names)
cell_names <- gsub("\\)", ".", cell_names)

outdir <- "/tscc/nfs/home/lebrusman/Gaulton_lab/code/pankbase_celltype_marker_genes/outputs_ND_only_250930/"

meta <- read.csv("merged_metadata.csv")

meta_singleCell <- meta[which(meta$treatments == "no_treatment" & meta$diabetes_status_description == "NonDiabetic"), ]


for (m in cells){
    message("Starting: ", m)
    #for pankbase
    m_name <- gsub("\\+", ".", m)
    m_name <- gsub("\\(", ".", m_name)
    m_name <- gsub("\\)", ".", m_name)

    meta_singleCell$samples <- gsub("-", ".", meta_singleCell$samples)
    keep_samps <- meta_singleCell$samples[which(meta_singleCell[,m_name]>=20)]

    interest <- read.table(paste0(dir, "cell_mtx/", m, '_persample_RNA_counts.tsv'), check.names = FALSE, header = TRUE, row.names = "gene")
    colnames(interest) <- gsub("-", ".", colnames(interest))
    other <- read.table(paste0(dir, "allbut_mtx/", m,'_persample_RNA_counts.tsv'), check.names = FALSE, header = TRUE, row.names = "gene")
    colnames(other) <- gsub("-", ".", colnames(other))
    
    interest <- interest[,colnames(interest) %in% keep_samps ]
    other <- other[,colnames(other) %in% keep_samps]

    colnames(interest) <- paste0('interest_', colnames(interest))
    colnames(other) <- paste0('other_', colnames(other))
    all <- cbind(interest, other)

    new_meta <- as.data.frame(colnames(all))
    colnames(new_meta) <- c('test')
    subtype <- str_extract(new_meta$test, "[^_]+")
    new_meta$subtype <- subtype

    new_meta$sample <- new_meta$test
    print(unique(new_meta$sample))

    new_meta$sample <- gsub('interest_', '', new_meta$sample)
    new_meta$sample <- gsub('other_', '', new_meta$sample)
    new_meta <- merge(new_meta, meta_singleCell, by.x='sample', by.y='samples')

   ff <- as.formula('~ subtype + sample') #Upate to your column with donor ID - maybe don't do this actually

   raw_counts <- all
   raw_counts <- raw_counts[,(colSums(raw_counts) > 0)] #Only test samples with any counts in the cell type
   meta_cell <- new_meta[which(new_meta$test %in% colnames(raw_counts)),]
   
   interest_raw_counts <- raw_counts[,grep('interest', colnames(raw_counts))]

    all_samps_cutoff <- floor(ncol(interest_raw_counts)/4)
    genes_to_keep <- list()
    for (i in 1:nrow(interest_raw_counts)) {
      if (sum(interest_raw_counts[i, ] >= 5) >= all_samps_cutoff) {
        genes_to_keep <- append(genes_to_keep, rownames(interest_raw_counts[i, ]))
      }
    }

    counts <- raw_counts
    counts <- counts[rownames(counts) %in% genes_to_keep,]

   rownames(meta_cell) <- meta_cell$test
   sample_order <- colnames(counts)
   meta_cell <- meta_cell[ sample_order,]

   my_design <- ff
   dds <- DESeqDataSetFromMatrix(round(counts), colData = meta_cell, design = my_design) #colData is where design columns are found
   dds <- estimateSizeFactors(dds)
   dds <- estimateDispersions(dds)
   dds2 <- dds
   print(dds2)

   ### Pairwise Wald test: conditon vs control  
   dds2 <- nbinomWaldTest(dds2)
   res <- results(dds2, contrast=c('subtype','interest','other'))
   res <- as.data.frame(res)
   res <- res[order(res$pvalue),]
   res <- cbind(gene = rownames(res), res)
   print(head(res))
   write.table(res, paste0(outdir, m,'_marker_genes_deseq.tsv'), sep='\t', quote=FALSE, row.names=FALSE,col.names=TRUE)

   #do fGSEA
   ### Run GSEA on each cell type and disease
    KEGG_react <- gmtPathways('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/reactome_kegg.gmt.txt')

    rpl <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rpl_file_gsea.csv', fill=TRUE, header=TRUE)
    rps <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rps_file_gsea.csv', fill=TRUE, header=TRUE)
    mtr <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/mts_file_gsea.csv', fill=TRUE, header=TRUE)

    ribo_proteins <- c(rpl$`Approved symbol`, rps$`Approved symbol`, mtr$`Approved symbol`)
    ribo_proteins <- ribo_proteins[which(ribo_proteins != 'Approved symbol')]
    
    
    res <- res[which(!rownames(res) %in% ribo_proteins),]
    res$pvalue <- ifelse(res$pvalue == 0, 1e-306, res$pvalue) 
    res$rank = (-log10(as.numeric(res$pvalue)))*res$log2FoldChange
    res = data.frame("SYMBOL" = rownames(res),
                     "stat" = res$rank)
    res = res[!grepl(pattern = "NA", x = res$SYMBOL),]
    res = res[!grepl(pattern = "MT-", x = res$SYMBOL),]
    ranks <- deframe(res)
    ranks = sort(ranks, decreasing=TRUE)
    deseq_df <- data.frame(gene = res$SYMBOL, GSEArank = res$stat)

    
    KEGG_react_fgseaRes <- fgseaMultilevel(pathways=KEGG_react,
                                stats=ranks,
                                eps =0.0,
                                minSize  = 0, 
                                maxSize  = 1000)#,nPermSimple = 10000)
    message("Number of total enriched terms ", m,": ", nrow(KEGG_react_fgseaRes))
    
    FDR_tresh = 0.10
    KEGG_react_fgseaRes.tresh = KEGG_react_fgseaRes[KEGG_react_fgseaRes$padj < FDR_tresh,]
    message("Number of significant terms: ", m,": ",nrow(KEGG_react_fgseaRes.tresh))
    ## Add categories
    res_signif <- KEGG_react_fgseaRes.tresh
    res_signif <- res_signif[order(res_signif$pval),]

    res_all <- KEGG_react_fgseaRes
    res_all <- res_all[order(res_all$pval),]
   
    fwrite(res_all, file=paste0(outdir, m, "_fGSEA_res_all.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)
    fwrite(res_signif, file=paste0(outdir, m, "_fGSEA_res_signif.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)
   
    message("Done with: ", m)

  } 


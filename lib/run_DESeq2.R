run_DESeq2 <- function(counts, meta, model, contrast, key) {
  DESeq2 <- lapply(contrast, function(cf) {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = model)
    
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    dds <- results(dds, filterFun = ihw, contrast = cf)
    dds <- as.data.frame(dds) |>
      rownames_to_column(var = "gene_id") |>
      left_join(key) |>
      relocate(symbol, .after = gene_id) |>
      as_tibble() |>
      arrange(padj)
    
  })
  
  # Generate comparison labels from contrasts  
  name <- unlist(contrast)
  name <- name[name!= "condition"]
  name <- split(as.character(name), 1:length(name) %% 2 == 0)
  name <- paste(name[[1]], name[[2]], sep = "_vs_")
  names(DESeq2) <- name
  return(DESeq2)
}
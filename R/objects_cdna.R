get.cdna.longest <- function(genes.df, genome, work.dir = NULL, root.name = NULL,
                             metas.include = c("cluster") # downstream of gr.find.cluster. remove unnecessary columns.
                             ) {
  # the only required column in genes.df is "gene" or "gene_name"
  if (is.character(genes.df)) {
    genes.df <- readRDS(genes.df)
  }
  
  if ("GRanges" %in% class(genes.df)) 
    genes.df <- utilsFanc::gr2df(genes.df)
  
  colnames(genes.df)[colnames(genes.df) == "gene_name"] <- "gene"
  if (!is.null(metas.include)) {
    genes.df <- genes.df[, colnames(genes.df) %in% c("gene", metas.include), drop = F]
  }
  
  utilsFanc::check.intersect("gene", "required column", colnames(genes.df), "colnames(genes.df)")
  
  if (is.null(work.dir)) {
    stop("work.dir must be specified")
  }
  
  cdnadb <- seqinr::read.fasta(paste0("~/genomes/", genome, "/ensembl/cdna/", genome, ".cdna.fa.gz"))
  headers <- unlist(seqinr::getAnnot(cdnadb))
  
  cdnadf <- data.frame(gene = stringr::str_extract(headers, "gene_symbol:[^ ]+ ") %>% gsub("^.+:| ", "", .),
                       seq = sapply(seqinr::getSequence(cdnadb), function(x) paste0(x, collapse = "")))
  

  dir.create(work.dir, showWarnings = F, recursive = T)
  
  cdnadf <- cdnadf %>% na.omit() %>% filter(gene %in% genes.df$gene)
  cdnadf$seqlength <- stringr::str_length(cdnadf$seq)
  
  cdnadf <- cdnadf %>% split(f = cdnadf$gene) %>% 
    lapply(function(df) {
      data.frame(gene = df$gene[1], seq = df$seq[which.max(df$seqlength)])
    }) %>% do.call(rbind, .)
  
  genes.df <- genes.df %>% dplyr::left_join(cdnadf)

  name <- "seqdf"
  if (is.null(root.name)) name <- paste0(root.name, "_", name)
  dir.create(work.dir, showWarnings = F, recursive = T)
  saveRDS(genes.df, paste0(work.dir, "/", name,".Rds"))
  
  invisible(genes.df)
}
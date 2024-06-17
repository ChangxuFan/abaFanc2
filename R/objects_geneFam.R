vec.find.cluster <- function(vec, max.gap = 5, min.length = 5,
                             bFlatten = F, sep.flatten = "_") {
  # vec = rep("a", 5) would return list(a = list(c(1,2,3,4,5)))
  # which means res$a[[1]]
  # if flatten, it would be res$a_1, if sep = "_".
  len <- length(vec)
  ir <- IRanges(start = 1:len, end = max.gap + (1:len), names = vec)
  irl <- split(ir, f = vec, drop = F)
  irl <- lapply(irl, function(ir) {
    ir.merge <- reduce(ir, drop.empty.ranges=TRUE, with.revmap=TRUE)
    revmap <- mcols(ir.merge)$revmap
    revmap.ori <- revmap
    for (i in 1:length(revmap.ori)) {
      revmap.ori[[i]] <- start(ir)[revmap.ori[[i]]]
    }
    mcols(ir.merge)$len <- sapply(revmap, length)
    mcols(ir.merge)$revmap.ori <- revmap.ori
    return(ir.merge)
  })
  irl.len <- lapply(irl, function(ir) {
    return(ir[mcols(ir)$len >= min.length])
  })
  
  res <- lapply(irl.len, function(ir) {
    l <- mcols(ir)$revmap.ori %>% as.list()
    return(l)
  }) %>% .[sapply(., length) > 0]
  if (bFlatten) {
    res <- lapply(names(res),function(name) {
      x <- res[[name]]
      names(x) <- paste0(name, sep.flatten, 1:length(x)) 
      return(x)
    }) %>% unlist(recursive = F)
  }
  return(res)
}

gr.find.cluster <- function(gr, name.col, cluster.col = "cluster",
                            max.gap = 5, min.length = 5, out.bed = NULL) {
  clusters <- vec.find.cluster(vec = mcols(gr)[, name.col],
                               max.gap = max.gap, min.length = min.length,
                               bFlatten = T)
  gr.clus <- lapply(names(clusters), function(cluster.name){
    gr.clus <- gr[(1:length(gr)) %in% clusters[[cluster.name]]]
    mcols(gr.clus)[, cluster.col] <- cluster.name
    return(gr.clus)
  }) %>% do.call(c, .)
  
  if (!is.null(out.bed)) {
    bed <- gr.clus %>% utilsFanc::gr2df() %>% 
      dplyr::mutate(cluster = !!as.name(cluster.col)) %>% 
      dplyr::select(chr, start, end, cluster, width, strand)
    collap <- bed %>%  dplyr::group_by(cluster) %>% 
      dplyr::summarise(chr = chr[1], start = min(start), end = max(end)) %>% 
      dplyr::ungroup() %>% dplyr::select(chr, start, end, cluster) %>% 
      as.data.frame()
    utilsFanc::write.zip.fanc(df = bed, out.file = out.bed, bed.shift = T)
    out.collap <- utilsFanc::insert.name.before.ext(out.bed, insert = "collap", delim = "_")
    utilsFanc::write.zip.fanc(df = collap, out.file = out.collap, bed.shift = T)
  }
  return(gr.clus)
}

gene.dist.cal <- function(seqdf, work.dir = NULL, root.name = NULL, first.n = NULL, threads = 1) {
  # first.n: for debug. Only work on the first n clusters
  utilsFanc::check.intersect(c("gene", "cluster", "seq"), "required columns", 
                             colnames(seqdf), "colnames(seqdf)")
  if (is.null(work.dir)) {
    work.dir <- tempdir()
  }
  dir.create(paste0(work.dir, "/mafft"), showWarnings = F, recursive = T)
  
  na.sum <- seqdf %>% dplyr::group_by(cluster) %>% dplyr::summarise(nonNA = sum(!is.na(seq))) %>% 
    dplyr::ungroup() %>% as.data.frame() %>% 
    dplyr::filter(nonNA > 1)
  
  if (!is.null(first.n)) {
    na.sum <- na.sum[1:min(nrow(na.sum), first.n),]
  }
  
  seqdf <- seqdf %>% dplyr::filter(cluster %in% na.sum$cluster)
  if (nrow(seqdf) < 1) {
    stop("nrow(seqdf) < 1")
  }

  mats <- seqdf %>% split(f = seqdf$cluster) %>% 
    utilsFanc::safelapply(function(seqdf) {
      # n.nonNA <- sum(!is.na(seqdf$seq))
      # if (n.nonNA < 2) {
      #   stop(paste0("Error at cluster ", cluster, ": at least 2 non-NA sequences needed"))
      # }
      seqdf <- seqdf[!is.na(seqdf$seq),]
      cluster <- seqdf$cluster[1]
      fa.in <- paste0(work.dir, "/mafft/", cluster, "_in.fa")
      fa.aligned <- paste0(work.dir, "/mafft/", cluster, "_aligned.fa")
      seqinr::write.fasta(sequences = as.list(seqdf$seq), names = seqdf$gene, file.out = fa.in, as.string = T)
      aln <- mafft.fanc(in.fa = fa.in, aligned.fa = fa.aligned)
      ape::dist.dna(ape::as.DNAbin(list(seq1 = c("-","-", "A", "T", "C", "-" ,"G"),
                                        seq2 = c("A", "A", "A", "C", "C", "-", "G"))), model = "raw")
      # gives 0.25. this means that gaps are not considered to begin with
      dist.mat <- ape::dist.dna(ape::as.DNAbin(aln), model = "raw", as.matrix = T)
      return(dist.mat)
    }, threads = threads)
  name <- "dist_mats"
  if (!is.null(root.name)) {
    name <- paste0(root.name, "_", name)
  }
  saveRDS(mats, paste0(work.dir, "/", name,".Rds"))
  return(mats)
}

gene.dist.summary <- function(dist.mats, work.dir = NULL, root.name = NULL) {
  # we first make sure there is no NA
  clusters <- names(dist.mats)
  na <- sapply(dist.mats, function(x) any(is.na(x)))
  if (sum(na) > 0) {
    warning(paste0("NA detected for ", sum(na), " clusters. First 5:\n",
                paste0(clusters[na][1:5], collapse = ", ")))
    dist.mats <- dist.mats[!na]
    clusters <- clusters[!na]
  }
  summary <- lapply(names(dist.mats), function(cluster) {
    mat <- dist.mats[[cluster]]
    stats <- list(cluster = cluster, n.genes = ncol(mat))
    
    diag(mat) <- 2
    flat <- as.vector(mat)
    names(flat) <- outer(rownames(mat), colnames(mat), function(x, y) return(paste0(x, ":", y))) %>% as.vector()
    flat <- flat[flat != 2]
    
    stats$min <- min(flat)
    stats$q1 <- quantile(flat, 0.25)
    stats$median <- median(flat)
    stats$q3 <- quantile(flat, 0.75)
    stats$max <- max(flat)
    
    stats$min.pair <- names(flat)[which.min(flat)]
    stats$q1.pair <- names(flat)[which.min(abs(flat - stats$q1))]
    stats$median.pair <- names(flat)[which.min(abs(flat - stats$median))]
    stats$q3.pair <- names(flat)[which.min(abs(flat - stats$q3))]
    stats$max.pair <- names(flat)[which.max(flat)]
    
    as.data.frame(stats)
  }) %>% do.call(rbind, .)
  rownames(summary) <- NULL
  if (!is.null(work.dir)) {
    dir.create(work.dir, showWarnings = F, recursive = T)
    name <- "dist_summary"
    if (!is.null(root.name)) 
      name <- paste0(root.name, "_", name)
    write.table(summary, paste0(work.dir,"/", name, ".tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
    
    pl <- lapply(c("min", "q1", "median", "q3", "max"), function(stat) {
      ggplot(summary, aes_string(x = stat)) + geom_histogram(bins = 30) + 
      ggtitle(stat)
    })
    scFanc::wrap.plots.fanc(pl, plot.out = paste0(work.dir, "/", name, "_histogram.png"))
  }
  invisible(summary)
}



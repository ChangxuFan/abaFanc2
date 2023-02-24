UGENE.bed.from.gtf <- function(loci.gr, gtf.type = "gene", gtf.col = "gene_name", out.file,
                               simplify.Ly49 = T) {
  # loci.gr: requires cols "regionID" and "gtf"
  if (length(gtf.col) != 1) {
    stop("length(gtf.col) != 1")
  }
  loci.gr$offset <- 1
  if (length(loci.gr) > 1)
    loci.gr$offset <- cumsum(c(1, width(loci.gr[1:(length(loci.gr)-1)])))
  df <- lapply(1:length(loci.gr), function(i) {
    gtf <- loci.gr$gtf[i] %>% rtracklayer::import()
    if ("type" %in% colnames(mcols(gtf))) {
      if (!gtf.type %in% gtf$type) {
        stop("!gtf.type %in% gtf$type")
      }
      gtf <- gtf[gtf$type %in% gtf.type]
    }
      
    mcols(gtf) <- mcols(gtf)[, gtf.col,drop = F]
    gtf <- subsetByOverlaps(gtf, loci.gr[i])
    start(gtf) <- start(gtf) - start(loci.gr)[i] + loci.gr$offset[i]
    end(gtf) <- end(gtf) - start(loci.gr)[i] + loci.gr$offset[i]
    
    df <- gtf %>% `names<-`(NULL) %>% as.data.frame()
    df$group <- loci.gr$regionID[i]
    df <- df[, c("seqnames", "start", "end", gtf.col, "group", "strand")]
    return(df)
  }) %>% do.call(rbind, .)
  if (simplify.Ly49) {
    df[, gtf.col] <- df[, gtf.col] %>% sub("^.+Ly49", "", .) %>% sub("^.+Gm44182.+$", "Gm", .)
  }
  utilsFanc::write.zip.fanc(df, out.file, bed.shift = F)
  invisible(df)
}


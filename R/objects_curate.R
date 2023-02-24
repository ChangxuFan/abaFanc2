rat.augment.exons <- function(query.df, ref.df, exons = 1:7,
                              buffer.left = 200, buffer.right = 200,
                              out.dir, threads = 10) {
  sub.df <- utilsFanc::safelapply(exons, function(i) {
    query <- query.df %>% filter(exon.number == i) %>% 
      dplyr::mutate(regionID = paste0("rn7.Ly49", gene),
                    genome = "rn7", id = forth, 
                    start = start - buffer.left, end = end + buffer.right) %>% 
      dplyr::select(chr, start, end, regionID, genome, strand, id)
    
    ref <- ref.df %>% filter(exon.number == i) %>% 
      dplyr::select(chr, start, end, strand, regionID, genome) %>% 
      dplyr::mutate(id = paste0("ref", i))
    aba.df <- rbind(ref, query)
    ext.dir <- paste0(out.dir, "/exon_", i, "_ext", buffer.left, "_", buffer.right)
    sub.dir <- paste0(out.dir, "/exon_", i, "_sub")
    abao1 <- aba.create(aba.df = aba.df, df.is.bed = F, 
                        work.dir = ext.dir)
    abao2 <- aba.subset(abao1, smart.cut.df = ref, new.work.dir = sub.dir)
    return(abao2@ori.df)
  }, threads = threads) %>%  Reduce(rbind, .)
  sub.df %>% split(., f = .$genome) %>% 
    lapply(function(df) {
      utilsFanc::write.zip.fanc(df[, c("chr", "start", "end", "id", "regionID", "strand")],
                                out.file = paste0(out.dir, "/all_sub_", df$genome[1], ".bed"), 
                                bed.shift = T)
    })
  return(sub.df)
}

augment.by.align <- function(query.df, ref.df, 
                             group.by = "exon", groups = NULL,
                              buffer.left = 200, buffer.right = 200,
                              out.dir, threads = 10, h = F) {
  cols.required <- c("chr", "start", "end", "strand", "genome", "regionID", group.by)
  
  if (h) {
    stop(paste0("cols required: ",
                paste0(cols.required, collapse = ", ")))
  }
  check.list <- list(query = query.df, ref = ref.df)
  lapply(names(check.list), function(name) {
    df <- check.list[[name]]
    not.found <- cols.required %>% .[! . %in% colnames(df)]
    if (length(not.found) > 0) {
      stop(paste0("required columns not found in ", name, ": ", 
                  paste0(not.found, collapse = ", ")))
    }
  })
  if (is.null(groups)) {
    groups <- intersect(unique(query.df[, group.by]), unique(ref.df[, group.by]))
  }
  if (length(groups) < 1) {
    stop("length(groups) < 1")
  }
  sub.df <- utilsFanc::safelapply(groups, function(group) {
    query <- query.df %>% filter(!!as.name(group.by) == group) %>% 
      dplyr::mutate(start = start - buffer.left, end = end + buffer.right) %>% 
      dplyr::select(chr, start, end, regionID, genome, strand, !!as.name(group.by))
    
    ref.cut <- ref.df %>% filter(!!as.name(group.by) == group) %>% 
      dplyr::select(chr, start, end, strand, regionID, genome, !!as.name(group.by))
    ref <- ref.cut %>% 
      dplyr::mutate(start = start - buffer.left, end = end + buffer.right)
    aba.df <- rbind(ref, query)
    ext.dir <- paste0(out.dir, "/", group, "_ext", buffer.left, "_", buffer.right)
    sub.dir <- paste0(out.dir, "/", group, "_sub")
    abao1 <- aba.create(aba.df = aba.df, df.is.bed = F, 
                        work.dir = ext.dir)
    abao2 <- aba.subset(abao1, smart.cut.df = ref.cut, new.work.dir = sub.dir)
    return(abao2@ori.df)
  }, threads = threads) %>%  Reduce(rbind, .)
  sub.df %>% split(., f = .$genome) %>% 
    lapply(function(df) {
      utilsFanc::write.zip.fanc(df[, c("chr", "start", "end", "regionID", "strand", group.by)],
                                out.file = paste0(out.dir, "/all_sub_", df$genome[1], ".bed"), 
                                bed.shift = T)
    })
  return(sub.df)
}

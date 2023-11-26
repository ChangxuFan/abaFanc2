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


paf.liftover <- function(bed, paf, min.l = 10000, lift.tempt.bed = NULL, lifted.bed = NULL, unlifted.bed = NULL,
                         return.lifted = F, return.dfs = F, 
                         paftools.js = PAFTOOLS.JS, add.params = "",
                         skip.paf = F) {
  root.name <- tools::file_path_sans_ext(bed)
  ext <- tools::file_ext(bed)
  if (is.null(lift.tempt.bed))
    lift.tempt.bed <- paste0(root.name, ".liftTempt.", ext)
  if (is.null(lifted.bed))
    lifted.bed <- paste0(root.name, ".lifted.", ext)
  if (is.null(unlifted.bed))
    unlifted.bed <- paste0(root.name, ".unlifted.", ext)
  cmd <- paste0(paftools.js, " liftover -l ", min.l, " ", add.params, " ", paf, " ", bed, " > ", lift.tempt.bed)
  print(cmd)
  if (skip.paf == F)
    system(cmd)
  
  lift.tempt <- read.table(lift.tempt.bed, as.is = T, header = F, comment.char = "") %>% 
    mutate(locus = sub("(_t[35])+$", "", V4)) %>% 
    select(V1, V2, V3, locus)
  if (return.lifted == T) {
    return(lift.tempt)
  }
  
  bed.df <- read.table(bed, as.is = T, header = F, comment.char= "") %>% 
    mutate(locus = paste0(V1, "_", V2, "_", V3))
  # locus <- paste0(bed.df$V1, "_", bed.df$V2, "_", bed.df$V3)
  

  lifted <- left_join(lift.tempt, bed.df[, 4:ncol(bed.df)])
  lifted <- lifted %>% .[,! colnames(.) %in% "locus"]
  
  unlifted <- bed.df[! bed.df$locus %in% lift.tempt$locus, ]
  
  trash <- utilsFanc::write.zip.fanc(lifted, out.file = lifted.bed)
  trash <- utilsFanc::write.zip.fanc(unlifted, out.file = unlifted.bed)
  res <- list(lifted.bed = lifted.bed, unlifted.bed = unlifted.bed, bed = bed, paf = paf)
  if (return.dfs == T) {
    res <- c(res, list(lifted.df = lifted, unlifted.df = unlifted, bed.df = bed.df, paf = paf))
  } 
  
  return(res)
}

paf.liftover.aotu <- function(main.paf.res = NULL, bed, paf, min.l = 10000, dist, size,
                              out.dir, 
                              paftools.js = PAFTOOLS.JS,
                              add.params = "",
                              debug = F, skip.paf = F) {
  stop("went through debugging but not fully tested yet!!")
  # ao means the Chinese character "ao". its opposite is "tu"
  dir.create(out.dir, showWarnings = F, recursive = T)
  root.name <- basename(bed)
  cmd <- paste0("ln -s ", normalizePath(bed, mustWork = T), " ", 
                normalizePath(out.dir, mustWork = T), "/", root.name)
  print(cmd); system(cmd)
  bed <- paste0(out.dir, "/", root.name)
  
  if (is.null(main.paf.res)) {
    main.paf.res <- paf.liftover(bed = bed, paf = paf, min.l = min.l,
                                 paftools.js = paftools.js, add.params = add.params,
                                 return.lifted = F, return.dfs = T,
                                 skip.paf = skip.paf)
    saveRDS(main.paf.res, paste0(out.dir, "/main.paf.res.Rds"))
  } else {
    if (is.character(main.paf.res)) {
      main.paf.res <- readRDS(main.paf.res)
    } else {
      saveRDS(main.paf.res, paste0(out.dir, "/main.paf.res.Rds"))
    }
  }
  utilsFanc::t.stat("done: lifting main bed")
  
  # main.paf.res$unlifted.bed <- "aotu/test_aotu/HG00741m.rm.simple.unlifted.bed"
  
  # generate the flanking region bed files
  # fu: up flank; fd: down flank
  fu.bed <- paste0(out.dir, "/", utilsFanc::insert.name.before.ext(name = root.name, insert = "fu", delim = "_"))
  fd.bed <- paste0(out.dir, "/", utilsFanc::insert.name.before.ext(name = root.name, insert = "fd", delim = "_"))
  
  cmd1 <- paste0("awk -F '\t' -v dist=",dist," -v size=",size,
                    " -v OFS='\t' '{$3 = $2 - dist; $2 = $3 - size;",
                    " if ($2 < 0) $2 = 0; if ($3 < 0) $3 = 0; print $0 }' ", 
                    main.paf.res$unlifted.bed, " > ", fu.bed)
  
  
  cmd2 <- paste0("awk -F '\t' -v dist=",dist," -v size=",size,
                " -v OFS='\t' '{$2 = $3 + dist; $3 = $2 + size;",
                " print $0 }' ", 
                main.paf.res$unlifted.bed, " > ", fd.bed)
  
  if (debug == T)
    threads <- 1
  else
    threads <- 2
  awk.res <- mclapply(list(cmd1, cmd2), function(cmd) {
    print(cmd)
    res <- system(cmd, intern = F)
    return(res)
  }, mc.cores = threads, mc.cleanup = T)
  if (any(awk.res != 0)) {
    stop("awk command failed")
  } else {
    utilsFanc::t.stat("done: taking flanking regions")
  }
  
  bed.list <- list(fu = fu.bed, fd = fd.bed)
  if (debug == T)
    threads <- 1
  else
    threads <- 2
  paf.res <- mclapply(names(bed.list), function(bed.name) {
    bed.file <- bed.list[[bed.name]]
    
    res <- paf.liftover(bed = bed.file, paf = paf, min.l = min.l,
                        paftools.js = paftools.js, add.params = add.params,
                        return.lifted = T,
                        skip.paf = skip.paf)
    return(res)
  }, mc.cores = threads)
  names(paf.res) <- names(bed.list)
  saveRDS(paf.res, paste0(out.dir, "/aotu.paf.res.Rds"))
  
  unlifted <- main.paf.res$unlifted.df %>%  
    mutate(start.pos = paste0(V1, "_", V2), end.pos = paste0(V1, "_", V3),
           fu = F, fd = F, aotu = F)
  
  fu.start <- parse.liftover.loci(paf.res$fu$locus) %>% 
    mutate(start = start + size + dist) %>% 
    mutate(start.pos = paste0(chr, "_", start)) %>% pull(start.pos)
  
  fd.end <- parse.liftover.loci(paf.res$fd$locus) %>% 
    mutate(end = end - size - dist) %>% 
    mutate(end.pos = paste0(chr, "_", end)) %>% pull(end.pos)
  
  unlifted$fu[unlifted$start.pos %in% fu.start] <- T
  unlifted$fd[unlifted$end.pos %in% fd.end] <- T
  
  unlifted$aotu[unlifted$fu == T & unlifted$fd == T] <- "aotu"
  unlifted$aotu[unlifted$fu == T & unlifted$fd == F] <- "fu.s"
  unlifted$aotu[unlifted$fu == F & unlifted$fd == T] <- "fd.s"
  unlifted$aotu[unlifted$fu == F & unlifted$fd == F] <- "flat"
  
  unlifted$aotu[(unlifted$V2 - dist - size) < 0 | (unlifted$V3 - dist - size) < 0] <- "trash"
  aotu <- unlifted %>% split(f = factor(unlifted$aotu, levels = c("aotu", "fu.s", "fd.s", "flat", "trash")))
  
  saveRDS(unlifted, paste0(out.dir, "/unlifted.Rds"))
  saveRDS(aotu, paste0(out.dir, "/aotu.Rds"))
  
  trash <- mclapply(names(aotu), function(x) {
    df <- aotu[[x]]
    trash <- utilsFanc::write.zip.fanc(df = df, out.file = paste0(out.dir, "/", x, ".bed"))
    return()
  })
  res <- list(aotu = aotu, unlifted.df = unlifted, 
              main.paf.res = paste0(out.dir, "/main.paf.res.Rds"),
              aotu.paf.res = paste0(out.dir, "/aotu.paf.res.Rds"))
  return(res)
}

parse.liftover.loci <- function(loci) {
  df <- data.frame(chr = gsub("_.+$", "", loci),
                   start = sub("(.+)_(.+)_(.+)", "\\2", loci) %>% as.numeric(),
                   end = sub("(.+)_(.+)_(.+)", "\\3", loci) %>% as.numeric())
  return(df)
}

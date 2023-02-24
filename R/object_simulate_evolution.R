simu.random.aln <- function(aln.file = NULL, ids, start = NULL, end = NULL, 
                            seq.list = NULL,
                            nuc.freq.list = NULL, length.list = NULL,
                            n.simu = 50,
                            # threads.pid = 1, 
                            threads.simu = 1, silent = F,
                            out.dir, root.name = NULL, mafft.params = MAFFT.AUTO,
                            ...) {
  # default: get nuc.freq and length from seq1 and seq2. But if you don't have sequences at hand, 
  #you can also just specify them
  # nuc.freq, if used, must be a named numeric vector of length 4 corresponding to ACGT 
  #and they must sum to 1! note:::: it's ACGT to conform to BioStrings' conventions! not ATCG!!
  library(Biostrings)
  if (silent) {
    do.plot <- do.save.rds <- write.pairwise.aln <- F
  } else {
    do.plot <- do.save.rds <- write.pairwise.aln <- T
  }
  if (is.null(root.name)) {
    if (!is.null(aln.file)) {
      root.name <- tools::file_path_sans_ext(basename(aln.file)) %>% paste0("_", start, "-", end)
    } else {
      root.name <- basename(out.dir)
    }
  }
  if (!is.null(aln.file)) {
    aln <- readDNAMultipleAlignment(aln.file)
    aln <- as.matrix(aln)
    if (is.null(start)) {
      start <- 1
    }
    if (is.null(end)) {
      end <- ncol(aln)
    }
  seq.list <- aln[ids, start:end] %>% t() %>% 
    as.data.frame() %>% `rownames<-`(NULL) %>% as.list()
    
  }
  # names(seqs) <- paste0(names(seqs), "@", start, "_", end)
  
  if (!is.null(seq.list)) {
    seq.list <- seq.list %>% lapply(function(x) {
      paste0(x, collapse = "") %>% gsub("-+", "", .) %>% return()
    }) 
    nuc.freq.list <- oligonucleotideFrequency(DNAStringSet(unlist(seq.list)), 1, as.prob = T) %>% 
      split(., f = 1:nrow(.))
    names(nuc.freq.list) <- names(seq.list)
    length.list <- seq.list %>% sapply(nchar)
  }
  trial.dir <- paste0(out.dir, "/", root.name, "_trials/")
  dir.create(trial.dir, showWarnings = F, recursive = T)
  pid.mats <- utilsFanc::safelapply(1:n.simu, function(simu.id) {
    bg.list <- lapply(1:length(seq.list), function(i) {
      freq <- nuc.freq.list[[i]]
      names(freq) <- c("A", "C", "G", "T")
      len <- length.list[[i]]
      bg <- lapply(names(freq), function(nuc) {
        f <- freq[nuc] * len %>% ceiling()
        return(rep(nuc, f))
      }) %>% Reduce(c, .)
      bg <- bg[1:len] %>% sample(size = len, replace = F)
      bg <- bg %>% paste0(collapse = "")
      # bg <- sample(c("A", "C", "G", "T"), size = len, replace = T,
      #              prob = freq) %>% paste0(collapse = "")
      return(bg)
    })
    names(bg.list) <- names(seq.list)#  %>% paste0("_rand", simu.id)
    aligned.fa <- paste0(trial.dir, "/", root.name, "_rand", simu.id, ".fa")
    trash <- mafft.fanc(in.fa = bg.list, aligned.fa = aligned.fa, mafft.params = mafft.params)
    pid.mat <- pairwise.dist.mat.ez(seq.mat = aligned.fa, model = "raw", use.pct.ident = T,
                                    pairwise.deletion = T, skip.plot = F)
    # pid.mat <- percent.identity.matrix(seq.vec = bg.list, out.dir = trial.dir,
    #                                    root.name = paste0(root.name, "_rand", simu.id),
    #                                    threads = threads.pid, 
    #                                    do.plot = do.plot, do.save.rds = do.save.rds, 
    #                                    write.pairwise.aln = write.pairwise.aln,
    #                                    ...)
    return(pid.mat)
  }, threads = threads.simu)
  fg.mat <- pairwise.dist.mat.ez(seq.mat = aln[ids, start:end], model = "raw", 
                                 use.pct.ident = T, pairwise.deletion = T, skip.plot = T)  
  # fg.mat <- percent.identity.matrix(seq.vec = seq.list, out.dir = trial.dir,
  #                                   root.name = paste0(root.name, "_fg"),
  #                                   threads = threads.pid, 
  #                                   do.plot = T, do.save.rds = T, 
  #                                   write.pairwise.aln = T, ...)
  # names(seq.list)
  pl <- utilsFanc::safelapply(names(seq.list), function(i) {
    utilsFanc::safelapply(names(seq.list), function(j) {
      bg.vec <- sapply(pid.mats, function(mat) return(mat[i, j]))
      fg <- fg.mat[i, j]
      e <- ecdf(bg.vec)
      pvalue <- round(1-e(fg), digits = 4)
      p <- ggplot(data.frame(x = bg.vec), aes(x = x)) + 
        geom_density() +
        geom_vline(xintercept = fg, color = "red") +
        xlim(0, 1) +
        theme(aspect.ratio = 1) +
        ggtitle(paste0(j, " ", i, " ", pvalue))
      return(p)
    }) %>% return()
  }) %>% Reduce(rbind, .)
  trash <- scFanc::wrap.plots.fanc(plot.list = pl, n.col = length(seq.list), sub.width = 3, sub.height = 3,
                          plot.out = paste0(out.dir, "/", root.name, "_distro.png"))
  return(list(fg = fg.mat, bg = pid.mats))
}

simulate.random.seq <- function(in.seq) {
  in.seq <- as.character(in.seq) %>% paste0(collapse = "")
  len <- nchar(in.seq)
  in.seq <- in.seq %>% gsub("-+", "", .)
  freq <- oligonucleotideFrequency(DNAStringSet(in.seq), 1, as.prob = T)[1, ] %>% as.vector()
  names(freq) <- c("A", "C", "G", "T")
  # A C G T
  bg <- lapply(names(freq), function(nuc) {
    f <- freq[nuc] * len %>% ceiling()
    return(rep(nuc, f))
  }) %>% Reduce(c, .)
  len2 <- min(len, length(bg))
  bg <- bg[1:len2] %>% sample(size = len2, replace = F)
  if (len > len2) {
    bg <- c(bg, rep("A", len - len2))
  }
  bg <- bg %>% paste0(collapse = "")
  return(bg)
}


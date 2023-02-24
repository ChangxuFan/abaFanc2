# gtf.add.intron <- function(gtf.gr, group.by = "transcript_id") {
#   warning("intron numbering changes with strand. Be careful...")
#   df <- gtf.gr %>% `names<-`(NULL) %>% as.data.frame()
#   intron.df <- df %>% dplyr::filter(type == "exon") %>% 
#     split(., f = .[, group.by]) %>% 
#     lapply(function(x) {
#       if (any(duplicated(x$exon_number))) {
#         stop("any(duplicated(x$exon_number))")
#       }
#       x <- x %>% dplyr::arrange(start)
#       y <- x
#       y$type <- "intron"
#       y <- y[1:(nrow(y)-1), ]
#       y$start <- x$end[-length(x$end)] + 1
#       y$end <- x$start[-1] - 1
#       # if (x$strand[1] == "-") {
#       #   browser()
#       #   y$exon_number <- (as.numeric(y$exon_number) - 1) %>% rev() %>% as.character()
#       #   # my desperate attempt to solve the strand issue that affects intron naming.
#       # }
#       return(y)
#     }) %>% Reduce(rbind, .)
#   all.gr <- rbind(df, intron.df) %>% dplyr::arrange(seqnames, start) %>% 
#     makeGRangesFromDataFrame(keep.extra.columns = T)
#   return(all.gr)
# }

gtf.add.intron <- function(gtf, group.by = "transcript_id",  out.file = NULL) {
  if (is.character(gtf)) {
    gtf <- rtracklayer::import(gtf)
  }
  gtf$exon_number <- as.numeric(gtf$exon_number)
  gtf$t.group <- mcols(gtf)[, group.by]
  
  df <- gtf %>% utilsFanc::gr2df()
  df <- df[df$type == "exon", ]
  df <- df %>% split(., f=.$t.group) %>% 
    lapply(function(df.sub) {
      dups <- df.sub$exon_number %>% .[duplicated(.)]
      if (length(dups) > 0) {
        stop(paste0("duplicated exon number for ", group.by, ": ", df.sub$t.group[1], "\n",
                    paste0(dups[1:min(5, length(dups))], collapse = ", ")))
      }
      intron.df <- lapply(df.sub$exon_number, function(e.this) {
        e.next <- e.this + 1
        if(! e.next %in% df.sub$exon_number) {
          return()
        }
        this.df <- df.sub %>% dplyr::filter(exon_number == e.this)
        next.df <- df.sub %>% dplyr::filter(exon_number == e.next)
        i.df <- this.df
        pos <- c(this.df$start, this.df$end, next.df$start, next.df$end) %>% sort()
        i.df$start <- pos[2]
        i.df$end <- pos[3]
        if (this.df$strand != next.df$strand) {
          stop("this.df$strand != next.df$strand")
        }
        i.df$type <- "intron"
        return(i.df)
      }) %>% do.call(rbind, .)
      df.sub <- rbind(df.sub, intron.df) 
      return(df.sub)
    }) %>% do.call(rbind, .)
  gtf <- df %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  names(gtf) <- NULL
  gtf$t.group <- NULL
  gtf <- gtf %>% gtf.fix.gene()
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    rtracklayer::export(gtf, out.file)
  }
  return(gtf)
}

gtf.filter.transcript <- function(gtf, group.genes.by = "gene_id", mode = "longest",
                                  out.file = NULL) {
  if (is.character(gtf)) {
    gtf <- rtracklayer::import(gtf)
  }
  if (!group.genes.by %in% colnames(mcols(gtf))) {
    stop("!group.genes.by %in% colnames(mcols(gtf))")
  }
  gtf$t.group <- mcols(gtf)[, group.genes.by]
  if (mode == "longest") {
    # gt: gene and transcript
    gt.df <- utilsFanc::gr2df(gtf)[, c("type","t.group", "transcript_id", "width")] %>% 
      dplyr::filter(type == "transcript") %>% 
      dplyr::group_by(t.group) %>% 
      dplyr::summarise(trans.pick = transcript_id[which.max(width)]) %>% 
      dplyr::ungroup()
    gtf <- gtf[is.na(gtf$transcript_id) | gtf$transcript_id %in% gt.df$trans.pick]
  } else {
    stop("only mode 'longest' has been developed")
  }
  
  gtf$t.group <- NULL
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    rtracklayer::export(gtf, out.file)
  }
  return(gtf)
}

gtf.fix.gene <- function(gtf, out.file = NULL, gene.id.col = "gene_name") {
  if (is.character(gtf)) {
    gtf <- rtracklayer::import(gtf)
  }
  if (! gene.id.col %in% colnames(mcols(gtf))) {
    stop("! gene.id.col %in% colnames(mcols(gtf))")
  }
  gtf$t.group <- mcols(gtf)[, gene.id.col]
  gtf <- gtf %>% sort()
  gtf <- gtf %>% split(f = factor(gtf$t.group, levels = unique(gtf$t.group))) %>%
    lapply(function(gr) {
      
      gr <- gr[gr$type != "gene"]
      add <- gr[1]
      add$exon_number <- NA
      add$transcript_id <- NA
      add$type <- "gene"
      start(add) <- min(start(gr))
      end(add) <- max(end(gr))
      gr <- c(add,gr)
      return(gr)
    }) %>% Reduce(c, .)
  gtf$t.group <- NULL
  
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    rtracklayer::export(gtf, out.file)
  }
  return(gtf)
}

gtf.align <- function(gtfs, genomes, entry.type = "exon", bFlipStrand = T,
                      split.by, splits = NULL, # eg: exon1, exon2...
                      region.by = "gene_name", regions = NULL,
                      out.dir, threads.aba = 1, npar = 1) {
  if (length(gtfs) != length(genomes)) {
    stop("length(gtfs) != length(genomes)")
  }
  gtfs <- lapply(1:length(gtfs), function(i) {
    gtf <- gtfs[[i]]
    if (is.character(gtf)) {
      gtf <- rtracklayer::import(gtf)
    }
    gtf$genome <- genomes[i]
    return(gtf)
  })
  gtf <- Reduce(c, gtfs)
  
  df <- gtf %>% utilsFanc::gr2df() %>% 
    dplyr::filter(type == entry.type)
  if(entry.type == "intron" && split.by == "exon_number") {
    df$intron_number <- df$exon_number
    split.by <- "intron_number"
  }
  
  if (!is.null(splits)) {
    df <- df[df[,split.by] %in% splits, ]
  }
  if (!is.null(regions)) {
    df <- df[df[,region.by] %in% regions, ]
  }
  df %>% split(., f = .[, split.by]) %>% 
    utilsFanc::safelapply(function(df.sub) {
      print(paste0("aligning ", split.by, ": ", df.sub[1, split.by]))
      aba.df <- df.sub %>% dplyr::mutate(regionID = !!as.name(region.by))
      
      if (bFlipStrand) {
        aba.df <- aba.df %>% dplyr::mutate(., strand = ifelse(.$strand == "+", "-", "+"))
      }
      aba.df <- aba.df[, c("chr", "start", "end", "regionID", "genome", "strand")]
      root <- paste0(split.by, "_", df.sub[1, split.by])
      if (entry.type %in% c("exon", "intron") && split.by %in% c("exon_number", "intron_number")) {
        root <- sub("exon_number_", "exon", root) %>% sub("intron_number_", "intron", .)
      }
      aba.create(aba.df, work.dir = paste0(out.dir, "/", root, "/"), df.is.bed = F, 
                 threads = threads.aba)
      return()
    }, threads = npar)
  return()
}


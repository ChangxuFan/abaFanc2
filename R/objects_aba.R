setClass(Class = "abao", slots = c(work.dir = "character", ori.df = "data.frame", ori.gr = "GRanges",
                                   ori.fa = "list", aln.fa = "list", map = "list", meta.data = "list"))

.DollarNames.abao <- function(x, pattern="") {
  grep(pattern, names(x@meta.data), value=TRUE)[1] %>% return()
}

"$.abao" <- function(x, i){
  val <- x@meta.data[[i]]
  return(val)
}

"$<-.abao" <- function(x, i, value){
  x@meta.data[[i]] <- value
  return(x)
}

"[.abao" <- function(x, i) {
  return(x@meta.data[i])
}
"[<-.abao" <- function(x, i, value) {
  x@meta.data[i] <- value
  return(x)
}

"[[<-.abao" <- function(x, i, value) {
  x@meta.data[[i]] <- value
  return(x)
}


aba.create <- function(aba.df, genome = NULL, fa = NULL, df.is.bed,
                       mask.lower.regex = NULL,
                       mafft = MAFFT.DEFAULT, mafft.params = MAFFT.AUTO,
                       do.prank = F, prank = PRANK,
                       do.phyml = F,
                       work.dir, threads = 1, fa.root.name = NULL) {
  if (is.null(fa.root.name))
    fa.root.name <- basename(work.dir) %>% paste0("_")
  # fields for aba.df: chr	start	end	regionID	genome/fa  strand.
  # strand seems okay to be missing.
  system(paste0("mkdir -p ", work.dir))
  work.dir <- normalizePath(work.dir, mustWork = T)
  if (is.character(aba.df)) {
    if (df.is.bed == T) {
      df.bed <- read.table(aba.df, sep = "\t", as.is = T)
      aba.df <- df.bed[, 1:4]
      colnames(aba.df) <- c("chr", "start", "end", "regionID")
      if (!is.null(df.bed$V6))
        aba.df$strand <- df.bed$V6
    } else {
      aba.df <- read.table(aba.df, sep = "\t", as.is = T, quote = "", header = T)
    }
  }

  aba.df$regionID <- gsub("[^A-Za-z0-9]", ".", aba.df$regionID)

  if (!is.null(genome))
    aba.df$genome <- genome
  if (!is.null(fa))
    aba.df$fa <- fa
  ori.fa <- get.fasta.bed.2(df = aba.df, genome = genome, fa = fa,
                            df.is.bed = df.is.bed, threads = threads)
  if (!is.null(mask.lower.regex)) {
    to.mask <- names(ori.fa) %>% .[grepl(paste0(mask.lower.regex, collapse = "|"), .)]
    if (length(to.mask) > 0) {
      ori.fa[to.mask] <- ori.fa[to.mask] %>% lapply(function(x) {
        x <- gsub("[a-z]", "N", x)
        return(x)
      })
    }
  }
  ori.fa.file <- paste0(work.dir, "/", fa.root.name, "in.fa")
  aln.fa <- mafft.fanc(in.fa = ori.fa, tempt.fa = ori.fa.file,
                      mafft = mafft, mafft.params = mafft.params,
                       aligned.fa = paste0(work.dir, "/", fa.root.name, "aligned.fa"))
  if (do.prank == T) {
    prank.root <- paste0(work.dir, "/prank/", basename(ori.fa.file) %>% sub("in.fa", "prank",.))
    trash <- prank.fanc(in.fa = ori.fa.file,
                        out.root.name = prank.root,
                        use.F = T, prank = prank)
  }

  map <- map.space(aligned.fa = aln.fa, threads = threads)
  abao <- new("abao", work.dir = work.dir, ori.df = aba.df,
              ori.gr = GenomicRanges::makeGRangesFromDataFrame(aba.df, keep.extra.columns = T),
              ori.fa = ori.fa, aln.fa = aln.fa, map = map,
              meta.data = list(mafft = mafft, mafft.params = mafft.params, fa.root.name = fa.root.name))

  fasta2phylip(in.fasta = paste0(abao@work.dir, "/", fa.root.name, "aligned.fa"))
  saveRDS(abao, paste0(abao@work.dir, "/abao.Rds"))
  aba.write.ori.df(abao)
  abao <- aba.add.consensus(abao)
  # abao <- aba.add.consensus(abao, shrink.all = T)
  if (do.phyml) {
    phyml.GTR.GIF(in.phy = paste0(abao@work.dir, "/", fa.root.name, "aligned.phy"))
  }
  invisible(abao)
}

map.space <- function(aligned.fa, threads =1) {
  # if (is.character(in.fa))
  #   in.fa <- read.fasta.fanc(in.fa)
  if (is.character(aligned.fa))
    aligned.fa <- read.fasta.fanc(aligned.fa)

  maps <- mclapply(aligned.fa, function(fa) {
    # browser()
    df <- data.frame(aln = 1:length(fa),
                     seq = fa)
    df2 <- df[df$seq != "-",]
    df2$ori <- 1:nrow(df2)
    df <- left_join(df, df2)
    return(df)
  }, mc.cores = threads, mc.cleanup = T)
  return(maps)
}

aba.subset <- function(abao, smart.cut.df = NULL, regionIDs = NULL, regionIDs.exclude = NULL,
                       ori.df.only = F,
                       new.work.dir, threads = 6,
                       check.grep.only = F,
                       mafft = MAFFT.DEFAULT, mafft.params = MAFFT.AUTO, ...) {
  # mode can be "union" or "intersect"
  # smart.cut.df is genome ordinates based. can be gr object.
  #the only metadata column used in regionID. buffer (which is present in a traditional smart.cut.df) is ignored
  # .cut.df: mapped to alignment space (the first nuc in aln space is 1.)
  warning("smart.cut.df only supports 1 region")
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }

  new.ori.df <- abao@ori.df
  if (!is.null(regionIDs)) {
    new.ori.df <- new.ori.df %>% dplyr::filter(grepl(paste0(regionIDs, collapse = "|"), regionID))
  }
  if (!is.null(regionIDs.exclude)) {
    new.ori.df <- new.ori.df %>% dplyr::filter(!grepl(paste0(regionIDs.exclude, collapse = "|"), regionID))
  }
  if (check.grep.only == T) {
    return(new.ori.df$regionID)
  }
  if (!is.null(smart.cut.df)) {
    .cut.df <- aba.map.2.consensus(abao = abao, df = smart.cut.df)
    new.ori.df <- new.ori.df %>% split(., f = factor(.$regionID, levels = .$regionID)) %>%
      lapply(function(df) {
        ori <- abao@map[[df$regionID]] %>% filter(aln >= .cut.df$start, aln <= .cut.df$end) %>%
          pull(ori) %>% .[!is.na(.)]
        if (length(ori) < 1)
          return()
        ori <- ori - 1
        if (!is.null(df$strand)) {
          if(df$strand == "-") {
            anchor <- df$end
            ori <- -1 * ori
          } else
            anchor <- df$start
        } else {
          anchor <- df$start
        }
        df$start <- anchor + min(ori)
        df$end <- anchor + max(ori)
        return(df)
      }) %>% Reduce(rbind, .)
  }
  if (ori.df.only) {
    return(new.ori.df)
  }
  abao <- aba.create(aba.df = new.ori.df, df.is.bed = F, work.dir = new.work.dir, threads = threads,
             mafft = mafft, mafft.params = mafft.params, ...)
  return(abao)
}


aba.subset.batch <- function(df, out.dir, threads.aba = 1, npar = 1) {
  # required columns: name	abao	ref.region	start	end	exclude	include
  cols <- c("name", "abao", "ref.region", "start", "end", "exclude", "include")
  if (is.character(df)) {
    df <- read.table(df, header = T)
  }
  not.found <- cols %>% .[!. %in%  colnames(df)]
  if (length(not.found) > 0) {
    stop(paste0("required columns not found: \n",
                paste0(not.found, collapse = ", "), "\n"))
  }
  df %>% split(f = 1:nrow(df)) %>%
    utilsFanc::safelapply(function(line) {
      print(paste0("processing line: ", line$name))
      smartcut <- readRDS(line$abao)@ori.df %>% .[.$regionID == line$ref.region,]
      if (nrow(smartcut) != 1) {
        stop("nrow(smartcut) != 1")
      }

      if (is.na(line$start)) {
        line$start <- 1
      }
      if (is.na(line$end)) {
        line$end <- smartcut$end - smartcut$start + 1
      }

      line$start <- max(line$start, 1)
      line$end <- min(line$end, smartcut$end)

      if (line$end < line$start) {
        stop("line$end < line$start")
      }

      if (!is.null(smartcut$strand) && smartcut$strand == "-") {
        warning("have not tested the negative strand scenario!")
        smartcut$end <- smartcut$end - line$start + 1
        smartcut$start <- smartcut$end - (line$end - line$start)
      } else {
        smartcut$start <- smartcut$start + line$start -1
        smartcut$end <- smartcut$start + line$end - line$start
      }

      out.dir <- paste0(out.dir, "/", line$name, "/")
      if (!is.na(line$include))
        regions.include <- unlist(strsplit(line$include, split = ", *"))
      else
        regions.include <- NULL

      if (!is.na(line$exclude))
        regions.exclude <- unlist(strsplit(line$exclude, split = ", *"))
      else
        regions.exclude <- NULL
      trash <- aba.subset(abao = line$abao, smart.cut.df = smartcut,
                          regionIDs = regions.include,
                          regionIDs.exclude = regions.exclude,
                          new.work.dir = out.dir, threads = threads.aba)
      return()
    }, threads = npar)
  return()
}



aba.subset.region.core <- function(abao, ori.df, .cut.df) {
  # format of .cut.df: chr (must be "aln"), start, end. other columns ignored.
  if (nrow(.cut.df) != 1)
    stop("nrow(.cut.df) != 1")
  new.ori.df <- ori.df %>% split(., f = factor(.$regionID, levels = .$regionID)) %>%
    lapply(function(df) {
      ori <- abao@map[[df$regionID]] %>% filter(aln >= .cut.df$start, aln <= .cut.df$end) %>%
        pull(ori) %>% .[!is.na(.)]
      if (length(ori) < 1)
        return()
      ori <- ori - 1
      if (!is.null(df$strand)) {
        if(df$strand == "-") {
          anchor <- df$end
          ori <- -1 * ori
        } else
          anchor <- df$start
      } else {
        anchor <- df$start
      }
      df$start <- anchor + min(ori)
      df$end <- anchor + max(ori)
      return(df)
    }) %>% Reduce(rbind, .)
  return(new.ori.df)
}

aba.consensus.2.bed <- function(abao, aln.df, out.dir, root.name = NULL, ignore.strand = F,
                                split.by = "genome", split.by.int.name = F) {
  # aln.df format: chr, start, end, int.name. chr must be "aln". other columns ignored
  # this is basically .cut.df. the only difference is that .cut.df allows only 1 line
  if (is.character(aln.df))
    read.table(aln.df, sep = "\t", header = T)
  aln.df$chr <- "aln"
  if (is.null(aln.df$int.name)) {
    aln.df$int.name <- "."
  }
  abao@ori.df %>% split(., f = .[, split.by]) %>%
    lapply(function(ori.df.sub) {
      split.name <- ori.df.sub[, split.by][1]
      bed <- aln.df %>% split(., f = 1:nrow(.)) %>%
        lapply(function(.cut.df) {
          bed <- aba.subset.region.core(abao = abao, ori.df = ori.df.sub, .cut.df = .cut.df)
          if (is.null(bed$strand))
            bed$strand <- "*"
          bed$name <- .cut.df$int.name
          if (ignore.strand == T) {
            return(bed[, c("chr", "start", "end", "name", "regionID")])
          }
          return(bed[, c("chr", "start", "end", "name", "regionID", "strand")])

        }) %>% Reduce(rbind, .)
      if (split.by.int.name == F) {
        out.file <- paste0(out.dir, "/", basename(abao@work.dir),
                           "_", root.name, "_", split.name, ".bed" )
        utilsFanc::write.zip.fanc(df = bed, out.file = out.file, bed.shift = T)
      } else {
        bed %>% split(., f = .$name) %>%
          lapply(function(bed.sub) {
            sub.name <- bed.sub$name[1]
            out.file <- paste0(out.dir, "/", basename(abao@work.dir),
                               "_", root.name, "_", split.name, "_", sub.name, ".bed" )
            utilsFanc::write.zip.fanc(df = bed.sub, out.file = out.file, bed.shift = T)
            return()
          })
      }
    })
}

aba.export.fasta <- function(abao, regions.df, mode, out.fa) {
  # mode can be "union" or "intersect"
  cut.df <- aba.map.2.consensus(abao = abao, df = regions.df)
  if (nrow(cut.df) > 1)
    stop("currently union or intersect are untested")
  if (nrow(cut.df) != 1)
    stop("internal check: nrow(cut.df) must be 1")
  if (cut.df$end <= cut.df$start) {
    stop("internal check: cut.df$end <= cut.df$start")
  }
  sub.aln <- abao@aln.fa %>% lapply(function(x) return(x[cut.df$start:cut.df$end]))
  dir.create(dirname(out.fa), showWarnings = F, recursive = T)
  seqinr::write.fasta(sequences = sub.aln, names = names(sub.aln), file.out = out.fa)
  return(out.fa)
}

aba.add.consensus <- function(abao, shrink.all = F,
                              ambiguityMap = "N", threshold = 0.5, regions.use = NULL,
                              consensus.name = NULL, consensus.doc = NULL, force = F) {
  if (is.null(consensus.name)) {
    if (identical(ambiguityMap, Biostrings::IUPAC_CODE_MAP))
      map.name <- "IUPAC"
    else
      map.name <- ambiguityMap
    consensus.name <- paste0("cons_", map.name, "_", threshold)
    if (shrink.all) {
      consensus.name <- paste0(consensus.name, "_A")
    }

  }
  if (is.null(consensus.doc)) {
    consensus.doc <- consensus.name
  }

  if (!is.null(abao@meta.data$consensus[[consensus.name]]) && force == F)
    stop(paste0(consensus.name, " already present"))
  seq.list <- abao@aln.fa
  if (!is.null(regions.use)) {
    regions.not.found <- regions.use %>% .[!.%in% names(seq.list)]
    if (length(regions.not.found) > 0)
      stop(paste0("some regions are not found: \n",
                  paste0(regions.not.found[1:5] %>% .[!is.na(.)], collapse = "\n")))
    seq.list <- seq.list[regions.use]
  }
  strings <- seq.list %>% sapply(function(x) {
    return(paste0(x, collapse = ""))
  })
  cons.mat <- Biostrings::consensusMatrix(x = strings)
  rownames(cons.mat) <- toupper(rownames(cons.mat))
  sub.mat <- Biostrings::nucleotideSubstitutionMatrix(4, -1, baseOnly = F)
  sub.mat <- cbind(rbind(sub.mat, "-"=-8), "-"=-8)
  # stolen from the msa examples. https://rdrr.io/bioc/msa/man/msaConservationScore-methods.html
  cons.score <- msa::msaConservationScore(cons.mat, sub.mat, gapVsGap=0)
  cons.string <- Biostrings::consensusString(strings, ambiguityMap = ambiguityMap,
                                             threshold = threshold)
  shrink.map <- data.frame(seq = cons.string %>% strsplit(split = "") %>% unlist(),
                           pos.full = 1:nchar(cons.string))

  if (shrink.all) {
    mat <- abao@aln.fa %>% as.data.frame() %>% t() %>% as.matrix()
    bmat <- !(mat %in% c("A", "T", "C", "G", "a", "t", "c", "g")) %>%
      matrix(nrow = nrow(mat))
    bGap <- colSums(bmat) > 0
    shrink.map <- shrink.map[!bGap, ]
  } else {
    shrink.map <- shrink.map %>% filter(seq != "-")
  }
  shrink.map$pos.shrink <- 1:nrow(shrink.map)
  cons.list <- list(seq = cons.string,
                    seq.shrink = cons.string %>% strsplit("") %>% unlist() %>%
                      .[shrink.map$pos.full] %>% paste0(collapse = ""),
                    cons.score = cons.score,
                    cons.score.shrink = cons.score[shrink.map$pos.full],
                    shrink.map = shrink.map,
                    name = consensus.name,
                    doc = consensus.doc)
  if (is.null(abao@meta.data$consensus))
    abao@meta.data$consensus <- list()
  abao@meta.data$consensus[[consensus.name]] <- cons.list

  out.dir <- paste0(abao@work.dir, "/cons/", consensus.name, "/")
  dir.create(out.dir, showWarnings = F, recursive = T)
  out.file <- paste0(out.dir, "/", abao@meta.data$fa.root.name, consensus.name, ".fa")
  seqinr::write.fasta(sequences =  cons.string, names = consensus.name, file.out = out.file)

  out.file <- paste0(out.dir, "/", abao@meta.data$fa.root.name, consensus.name, "_shrink.fa")
  seqinr::write.fasta(sequences =  cons.list$seq.shrink, names = consensus.name, file.out = out.file)

  lapply(c(T, F), function(bOmit.cons) {
    lapply(c(T, F), function(bShrink) {
      omit <- ifelse(bOmit.cons, "",  "_wCons")
      shrink <- ifelse(bShrink, "_shrink", "")
      out.file <- paste0(out.dir, "/", abao@meta.data$fa.root.name, consensus.name,
                         shrink, omit, ".fa")
      aba.write.consensus.with.aln.2(abao = abao, cons.name = consensus.name,
                                     omit.cons = bOmit.cons, shrink = bShrink,
                                     out.file = out.file)
    })
  })

  # write conservation score as bdg:
  lapply(c("cons.score", "cons.score.shrink"), function(x) {
    vec <- cons.list[[x]]
    l <- length(vec)
    if (grepl("shrink", x))
      chr <- consensus.name %>% paste0("_shrink")
    else
      chr <- consensus.name
    df <- data.frame(chr = chr, start = (1:l) - 1, end = 1:l, score = vec)
    trash <- utilsFanc::write.zip.fanc(df = df,
                                       out.file = paste0(out.dir, "/", abao@meta.data$fa.root.name,
                                                         consensus.name, "_", x,".bdg"),
                                       bed.shift = F)
    return()
  })
  return(abao)
}

aba.write.consensus.with.aln <- function(abao, cons.name,  smart.cut.df = NULL, shrink = F, out.file = NULL) {
  if (shrink == T)
    cons.string <- abao@meta.data$consensus[[cons.name]]$seq.shrink
  else
    cons.string <- abao@meta.data$consensus[[cons.name]]$seq
  if (is.null(cons.string)) {
    stop(paste0("consensus named ", cons.name, " is not found"))
  }

  strings  <- abao@aln.fa %>% sapply(function(x) {
    if (shrink == T)
      x <- x[abao@meta.data$consensus[[cons.name]]$shrink.map$pos.full]
    return(paste0(x, collapse = ""))
  })

  strings <- c(cons.string, strings)
  names(strings)[1] <- cons.name
  if (is.null(out.file)) {
    out.file <- paste0(abao@work.dir, "/", basename(abao@work.dir), "_", cons.name,"_", shrink,"_wAln.fa")
  }
  dir.create(dirname(out.file), showWarnings = F, recursive = T)

  seqinr::write.fasta(sequences = as.list(strings), names = names(strings),
                      file.out = out.file, as.string = F)
  return()
}

aba.write.consensus.with.aln.2 <- function(abao, cons.name, aln.list = NULL,
                                           regionIDs.include = NULL, regionIDs.exclude = NULL,
                                           shrink = F,  smart.cut.df = NULL, omit.cons = F,
                                           max.gap.frac = 1,
                                           out.file = NULL) {

  if (is.null(aln.list)) {
    if (is.character(abao)) {
      abao <- readRDS(abao)
    }
    aln.df <- abao@aln.fa %>% as.data.frame()
    names(aln.df) <- names(abao@aln.fa)
    cons.string <- abao@meta.data$consensus[[cons.name]]$seq
    if (is.null(cons.string)) {
      stop(paste0("consensus named ", cons.name, " is not found"))
    }
    aln.df[, cons.name] <- cons.string %>% strsplit(split = "") %>% unlist()
    aln.df$pos.full <- 1:nrow(aln.df)
    aln.df <- left_join(aln.df, abao@meta.data$consensus[[cons.name]]$shrink.map)
    if (!is.null(smart.cut.df))
      aln.df <- aba.smart.cut(abao = abao, track.bp.df = aln.df, smart.cut.df = smart.cut.df,
                              location.col = "pos.full", n.breaks = NULL)[[1]]
    if (shrink == T)
      aln.df <- aln.df[!is.na(aln.df$pos.shrink),]

    aln.list <- aln.df[, c(cons.name, names(abao@aln.fa))] %>% as.list()
    if (omit.cons) {
      aln.list <- aln.list[-1]
    }
  }

  if (max.gap.frac < 1) {
    bPass <- sapply(aln.list, function(x) {
      t <- table(x)
      if (is.na(t["-"])) {
        n.gap <- 0
      } else {
        n.gap <- t["-"]
      }
      res <- (n.gap/sum(t)) < max.gap.frac
      return(res)
    }) %>% `names<-`(NULL)
    aln.list <- aln.list[bPass]
  }
  if (!is.null(regionIDs.include)) {
    aln.list <- aln.list[grepl(paste0(regionIDs.include, collapse = "|"), names(aln.list))]
  }
  if (!is.null(regionIDs.exclude)) {
    aln.list <- aln.list[!grepl(paste0(regionIDs.exclude, collapse = "|"), names(aln.list))]
  }

  if (is.null(out.file)) {
    out.file <- paste0(abao@work.dir, "/", basename(abao@work.dir), "_", cons.name,"_", shrink,"_wAln.fa")
  }
  dir.create(dirname(out.file), showWarnings = F, recursive = T)

  seqinr::write.fasta(sequences = aln.list, names = names(aln.list),
                      file.out = out.file, as.string = F)
  return()

}

aba.write.consensus <- function(abao, cons.name = NULL, shrink, out.dir = NULL, root.name = NULL) {
  if (is.null(cons.name)) {
    if (length(abao@meta.data$consensus) != 0) {
      stop("cons.name is null, and length(abao@meta.data$consensus) != 0")
    }
    cons.name <- names(abao@meta.data$consensus)
  }
  if (shrink == T) {
    slot <- "seq.shrink"
    suffix <- "_shrink"
  } else {
    slot <- "seq"
    suffix <- ""
  }

  if (is.null(out.dir)) {
    out.dir <- abao@work.dir
  }

  system(paste0("mkdir -p ", out.dir))

  if (is.null(root.name)) {
    root.name <- basename(abao@work.dir)
  }
  outfile <- paste0(out.dir, "/", root.name,"_", cons.name, suffix, ".fa")

  seqinr::write.fasta(sequences = abao@meta.data$consensus[[cons.name]][[slot]],
                      as.string = T, nbchar = 80, names = cons.name,
                      file.out = paste0(outfile))
  return()
}

bp.df.gen <- function(data.gr, region.df=NULL) {
  if (is.null(region.df))
    region.df <- data.gr
  if (!is.data.frame(region.df)) {
    region.df <- region.df %>% `names<-`(NULL) %>% as.data.frame()
    region.df$chr <- region.df$seqnames
  }

  start <- region.df$start %>% min()
  end <- region.df$end %>% max()

  if (length(seqnames(data.gr) %>% unique()) > 1)
    stop("only one chromosome is allowed")
  if (length(region.df$chr %>% unique()) > 1)
    stop("only one chromosome is allowed")
  if (length(seqnames(data.gr)) > 0) {
    if (as.character(region.df$chr[1]) != as.character(seqnames(data.gr)[1]))
      stop("chromosome must match")
  }
  chr <- as.character(region.df$chr[1])
  bp.gr <- data.frame(chr = chr, start = start:end, end = start:end) %>%
    GenomicRanges::makeGRangesFromDataFrame()
  j.gr <- bp.gr %>% plyranges::join_overlap_left(data.gr)
  return(j.gr)
}

aba.add.track <- function(abao, track.df = NULL, bTrack.df.regionID.regex = F,
                          track.file, track.name, track.regionID.regex = NULL, prepend.regionID = F, threads = 1,
                          track.type, smooth.win.1side = NULL) {
  # required fields: regionID	track.name	track.type	track.file
  if (is.null(track.df)) {
    if (is.null(track.regionID.regex)) {
      track.regionID <- abao@ori.df$regionID
    } else {
      track.regionID <- abao@ori.df$regionID %>% .[grepl(paste0(track.regionID.regex, collapse = "|"), .)]
    }
    if (is.null(track.name)) {
      track.name <- track.type
    }
    track.df <- data.frame(regionID = track.regionID, track.name = track.name,
                           track.type = track.type, track.file = track.file)
    prepend.regionID <- T
  }
  if (is.character(track.df))
    track.df <- read.table(track.df, header = T, as.is = T, sep = "\t", quote="")
  # browser()
  if (bTrack.df.regionID.regex) {
    track.df <- track.df %>% split(., f = 1:nrow(.)) %>%
      lapply(function(track.df) {
        track.df <- cbind(abao@ori.df$regionID %>% .[grepl(track.df$regionID[1], .)], track.df[, colnames(track.df) != "regionID"])
        colnames(track.df)[1] <- "regionID"
        return(track.df)
      }) %>% Reduce(rbind, .)
  }
  track.df$regionID <- gsub("[^A-Za-z0-9]", ".", track.df$regionID)
  track.df$track.name <- gsub("[^A-Za-z0-9]", ".", track.df$track.name)
  track.df <- track.df %>% filter(regionID %in% abao@ori.df$regionID)

  if (prepend.regionID == T) {
    track.df$track.name <- paste0(track.df$regionID, "..", track.df$track.name)
  }
  abao$track.df <- track.df
  if (sum(duplicated(track.df$name)) > 0) {
    stop("track name must be unique")
  }
  track.bp <- track.df %>% split(., f = factor(.$track.name, levels = .$track.name)) %>%
    mclapply(function(track) {
      print(track$track.name)
      region <- abao@ori.gr %>% plyranges::filter(regionID == track$regionID)
      # browser()
      if (track.type == "bdg") {
        # if (region$regionID == "129.Ly49e") {
        #   browser()
        # }
        if (grepl(".bdg(.gz)*$", track$track.file[1])) {
          gr <- rtracklayer::import.bedGraph(con = track$track.file, which = region)
        } else {
          gr <- rtracklayer::import(con = track$track.file, which = region)
        }

        if (!is.null(smooth.win.1side)) {
          gr <- utilsFanc::gr.expand.smooth(gr, smooth.win.1side = smooth.win.1side,
                                            drop.seqlevel = T, drop.zero = T)
        }
        gr.bp <- bp.df.gen(data.gr = gr, region.df = region)
        ori <- abao@ori.df %>% filter(regionID == track$regionID)
        if ("strand" %in% colnames(ori) && ori$strand == "-") {
          gr.bp$score <- rev(gr.bp$score)
        }
        if (is.null(gr.bp$score))
          gr.bp$score <- NA
        gr.bp$score[is.na(gr.bp$score)] <- 0
      } else if (track.type %in% c("bed", "rmsk")) {
        gr <- import.bed.fanc(bed = track$track.file, return.gr = T)
        gr.bp <- bp.df.gen(data.gr = gr, region.df = region)
        if (is.null(gr.bp$forth))
          gr.bp$forth <- NA
        if (is.null(gr.bp$fifth))
          gr.bp$fifth <- NA
      } else {
        stop("only bedgraph, bed, and rmsk are supported right now")
      }
      return(gr.bp)
    }, mc.cleanup = T, mc.cores = threads)

  abao[[paste0("track.bp.", track.type)]] <- track.bp
  track.bp.df <- mclapply(seq_along(track.bp), function(i) {
    # track <- track.bp$Klra8..CAGE
    track.name.i <- names(track.bp)[i]
    track <- track.bp[[i]] %>% `names<-`(NULL) %>% mcols() %>% as.data.frame()
    track <- track %>% mutate(ori = 1:nrow(track), track.name = track.name.i) # %>% select(ori, score)

    regionID <- abao$track.df %>% filter(track.name == track.name.i) %>% pull(regionID)
    map <- abao@map[[regionID]]
    track.mapped <- left_join(track, map) # %>% pull(score)
    if (track.type == "bdg") {
      track.mapped <- track.mapped %>% filter(!is.na(score))
    } else if (track.type %in% c("bed", "rmsk")) {
      track.mapped <- track.mapped %>% filter(!is.na(forth))
    }
    return(track.mapped)
  }, mc.cores = threads, mc.cleanup = T) %>% Reduce(rbind, .)
  # names(track.bp.df) <- names(track.bp)
  # track.bp.df <- as.data.frame(track.bp.df)
  #browser()
  abao[[paste0("track.bp.df.",track.type)]]  <- track.bp.df
  return(abao)
}

get.dendro.order <- function(df, axis) {
  if (axis == 2)
    df <- t(df)
  dend <- hclust(dist(df, diag = T))
  order <- dend$labels[dend$order]
  return(order)
}

aba.write.ori.df <- function(abao, as.bed = T, split.by = "genome") {
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  if (as.bed == T) {
    gr.list <- abao@ori.gr %>% split(., f = mcols(.)[, split.by])
    trash <- lapply(names(gr.list), function(split.name) {
      gr <- gr.list[[split.name]]
      out.file <- paste0(abao@work.dir, "/ori.df.bed/", basename(abao@work.dir), "_ori_gr_", basename(split.name), ".bed")
      utilsFanc::write.zip.fanc(df = gr, bed.shift = T, out.file = out.file)
      return()
    })
    return()
  } else {
    stop("not developed. currently only output bed files")
  }

}

aba.plot.hm <- function(abao, track.type, smart.cut.df = NULL, tracks.include = NULL,
                        use.mafft.order = F, use.order = NULL,
                        broadcast = F, fill.gap = F, bed.same.color = F,
                        abs = F, scale.row = F, normalize.row = T,normalize.to.max = F,
                        cluster.tracks = F,
                        remove.zero.tracks = F,
                        x.lim = NULL, project.x.to = NULL,
                        add.grid = T, grid.line.size = 0.1,
                        plot.out = NULL, height=5, width=20, text.size = NULL,
                        add.nucleotide = F, gap.only = F,
                        atcg.size = NULL, sub.align.out = NULL) {
  # arguments specific for bdg tracks: abs, scale.row, cluster.tracks, remove.zero.tracks
  track.bp.df <- abao@meta.data[[paste0("track.bp.df.", track.type)]]
  if (!is.null(x.lim)) {
    track.bp.df <- track.bp.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
  }

  if (!is.null(smart.cut.df)) {
    cut.res <- aba.smart.cut(abao = abao, track.bp.df = track.bp.df, smart.cut.df = smart.cut.df,
                             project.breaks.to = project.x.to)
    track.bp.df <- cut.res[[1]]
  }

  if (!is.null(tracks.include)) {
    track.bp.df <- track.bp.df %>% filter(track.name %in% tracks.include)
  }

  tracks.all <- track.bp.df %>% pull(track.name) %>% unique()
  tracks <- tracks.all
  # track.bp.df <- track.bp.df %>% filter(track.name %in% tracks)
  if (!is.null(use.order)) {
    tracks <- use.order %>% .[. %in% tracks]
    if (length(tracks) < 1) {
      stop("length(tracks) < 1")
    }
  } else if (use.mafft.order == T) {
    tracks <- data.frame(regionID = names(abao@map)) %>% left_join(abao$track.df) %>% pull(track.name) %>%
      .[!is.na(.)]
  }
  if (track.type == "bdg") {
    if (remove.zero.tracks == T) {
      track.sum <- track.bp.df %>% group_by(track.name) %>% filter(!is.na(score)) %>% summarise(sum = sum(score))
      tracks.non.zero <- track.sum %>% filter(sum != 0) %>% pull(track.name) %>% unique()
      tracks <- tracks[tracks %in% tracks.non.zero]
    }
    track.bp.df <- track.bp.df %>% filter(track.name %in% tracks)
    if (abs == T) {
      track.bp.df$score <- abs(track.bp.df$score)
    }

    if (cluster.tracks == T && use.mafft.order == F && is.null(use.order)) {
      non.all.zero <- track.bp.df %>% na.omit() %>% group_by(aln) %>%
        summarise(nz = sum(score != 0) > 0) %>% filter(nz == T) %>% pull(aln)

      track.scale.center <- track.bp.df %>% filter(aln %in% non.all.zero) %>%
        group_by(track.name) %>% mutate(score = scale(score, center = T)) %>%
        ungroup() %>% as.data.frame()
      mat <- reshape2::acast(track.scale.center, formula = aln ~ track.name, value.var = "score")
      tracks <- get.dendro.order(as.data.frame(mat), axis = 2)
    }

    if (scale.row == T) {
      track.bp.df <- track.bp.df %>% group_by(track.name) %>% mutate(score = scale(score, center = F)) %>%
        ungroup() %>% as.data.frame()
    } else if (normalize.row == T) {
      if (normalize.to.max) {
        track.bp.df <- track.bp.df %>% group_by(track.name) %>% mutate(score = score/(max(score) + 1)) %>%
          ungroup() %>% as.data.frame()
      } else {
        track.bp.df <- track.bp.df %>% group_by(track.name) %>% mutate(score = score/(sum(score) + 1)) %>%
          ungroup() %>% as.data.frame()
      }

    }

    if (!is.null(smart.cut.df)) {
      track.bp.df$aln <- factor(track.bp.df$aln, levels = cut.res$levels)
    }
    col.fun <- circlize::colorRamp2(c(0, 1), c("white", "red4"))
    p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = score)) +
      geom_tile() +
      # scale_fill_gradient(low="white", high="red4") +
      #scale_fill_gradientn(colours = c("white", "red1", "red2", "red4"),
      #                     values = c(0, 0.1, 0.8, 1.0)) +
      scale_fill_gradientn(colours = c("white", col.fun(0.15), "red4"),
                           values = c(0, 0.05, 1.0)) +
      theme_classic() +
      theme(panel.background = element_rect(fill = "grey88",
                                            colour = "grey88")) +
      theme(legend.position = "bottom")
    # browser()
  }

  if (track.type %in% c("bed", "rmsk")) {
    if (fill.gap == T) {
      #browser()
      track.bp.df <- track.bp.df %>% split(., f= factor(.$track.name, levels = unique(.$track.name))) %>%
        lapply(function(df) {
          regionID <- aba.get.regionID(abao, df$track.name[1])
          feature.bound <- df %>% group_by(forth) %>% summarise(start = min(aln), end = max(aln)) %>%
            ungroup() %>% as.data.frame()
          j <- left_join(abao@map[[regionID]], df)
          filled <- feature.bound %>% split(., f=1:nrow(.)) %>%
            lapply(function(feature) {
              filled.sub <- j %>% filter(aln >= feature$start, aln <= feature$end)
              filled.sub$forth <- feature$forth
              filled.sub$track.name <- df$track.name[1]
              return(filled.sub)
            }) %>% Reduce(rbind, .)
          return(filled)
        }) %>% Reduce(rbind,.)
    }

    if (broadcast == T) {
      track.bp.df$track.name <- "consensus"
      track.bp.df <- track.bp.df %>% unique()
      tracks <- "consensus"
    }

    if (!is.null(smart.cut.df)) {
      track.bp.df$aln <- factor(track.bp.df$aln, levels = cut.res$levels)
    }

    if (track.type == "bed") {
      if (bed.same.color) {
        p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks))) +
          geom_tile(fill = "#0000b2") +
          theme_classic()
      } else {
        p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = forth)) +
          geom_tile() +
          theme_classic()

      }

    }

    if (track.type == "rmsk") {
      track.bp.df$fifth <- sub("/.+$", "", track.bp.df$fifth)
      p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = fifth)) +
        geom_tile() +
        scale_fill_manual(values = c("darkorange", "grey20", "seagreen4", "grey30", "grey40", "red", "grey50", "blue"),
                          breaks = c("LINE", "Low_complexity", "LTR", "Satellite", "Simple_repeat", "SINE", "Unknown", "DNA"))
    }
  }

  if (!is.numeric(track.bp.df$aln)) {
    vlines <- cut.res$cut.df.buffer$start

    p <- p + scale_x_discrete(drop = F, breaks = factor(cut.res$breaks[-1], # %>% c(vlines),
                                                        levels = cut.res$levels),
                              labels = cut.res$breaks.proj[-1]) +
      geom_vline(xintercept = factor(vlines, levels = cut.res$levels)[-1], linetype = "dashed")

  }
  if (!is.null(text.size)) {
    p <- p + theme(text = element_text(size = text.size))
  }
  if (add.grid) {
    n.hlines <- track.bp.df$track.name %>% unique() %>% length()
    hlines <- 1:(n.hlines-1) + 0.5
    p <- p + geom_hline(yintercept = hlines, linetype = "dashed", size = grid.line.size)
  }

  if (add.nucleotide == T) {
    p <- aba.plot.aln(abao, tracks = tracks, p.in = p, x.lim = x.lim, gap.only = gap.only,
                      smart.cut.df = smart.cut.df, atcg.size = atcg.size, sub.align.out = sub.align.out)
  }
  p <- p + theme(axis.line = element_blank(), axis.ticks.y = element_blank())

  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = F)
    ggsave(plot.out, p, width = width, height = height, units = "in", dpi = 100, limitsize = F)
  }
  return(p)

}

aba.write.track <- function(abao, smart.cut.df = NULL, track.name = NULL, track.type,
                            normalize = F, scale.to = 1000,
                            add.constant = NULL, add.per.site = T,
                            tracks.include = NULL,
                            cluster.tracks = F, cluster.smart.cut = NULL,
                            track.order = NULL, region.as.sort.key = T, sort.missing.method = "add", order.only = F,
                            broadcast = F, fill.gap = T,
                            consensus.name, shrink = T,
                            write.consensus = T, consensus.add.N = NULL,
                            x.lim = NULL,
                            return.matrix = F, out.dir = NULL, root.name = NULL,
                            samtools = SAMTOOLS, bedtools = BEDTOOLS,
                            push.2.igv = T) {
  # note: track.name used to be referred to as track.type
  dir.create(out.dir, recursive = T, showWarnings = F)
  if (track.type == "bdg")
    col <- "score"
  else if (track.type %in% c("rmsk", "bed"))
    col <- "forth"
  else
    stop("track type not supported")

  if (is.null(root.name))
    root.name <- abao@meta.data$fa.root.name
  if (is.null(track.name))
    track.name <- track.type
  track.bp.df <- abao@meta.data[[paste0("track.bp.df.", track.name)]]
  if (is.null(track.bp.df))
    stop("is.null(track.bp.df)")

  if (!is.null(x.lim)) {
    track.bp.df <- track.bp.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
  }

  if (!is.null(smart.cut.df)) {
    cut.res <- aba.smart.cut(abao = abao, track.bp.df = track.bp.df, smart.cut.df = smart.cut.df)
    track.bp.df <- cut.res[[1]]
  }

  if (shrink == T) {
    shrink.map <- abao@meta.data$consensus[[consensus.name]]$shrink.map
    if (is.null(shrink.map))
      stop("is.null(shrink.map)")
    names(shrink.map) <- c("cons.seq", "aln", "pos.shrink")
    if (!is.null(smart.cut.df)) {
      shrink.map <- cut.res$cut.df.buffer %>% split(., f = 1:nrow(.)) %>%
        lapply(function(cut) {
          df <- shrink.map %>% filter(aln > cut$start, aln < cut$end)
          return(df)
        }) %>% Reduce(rbind, .)
    }
    shrink.map$pos.out <- 1:nrow(shrink.map)
    ref.fa <- paste0(out.dir, "/", root.name, "shrink.fa")
    seq.name <- paste0(consensus.name, "_shrink")
    if (write.consensus == T && !is.null(out.dir)) {
      seq.tmp <- shrink.map$cons.seq
      if (!is.null(consensus.add.N)) {
        seq.tmp <- c(seq.tmp, rep("N", consensus.add.N))
      }
      seqinr::write.fasta(sequences = list(seq.tmp),
                          names = seq.name,
                          file.out = ref.fa)

      system(paste0(samtools, " faidx ", ref.fa))
      utilsFanc::bash2ftp(filename = ref.fa)
    }
    ref.length <- length(unique(shrink.map$pos.out))
    track.bp.df <- inner_join(shrink.map, track.bp.df)
    if (track.type == "bdg" && !return.matrix)
      track.bp.df <- track.bp.df %>% filter(score != 0)

    if (track.type %in% c("bed", "rmsk")) {
      track.bp.df <- track.bp.df %>% filter(!is.na(forth))
      if (fill.gap == T) {
        track.bp.df <- aba.track.fill.gap.2(abao = abao, track.bp.df = track.bp.df, fill.col = "pos.out") %>%
          filter(!is.na(pos.out))
      }

      if (broadcast == T) {
        track.name <- track.bp.df$track.name %>% unique()
        if (length(track.name) > 1) {
          track.name <- track.name %>% sub("^.+\\.\\.", "", .) %>% unique()
          if (length(track.name) > 1)
            stop("length(track.name) > 1")
        }
        track.bp.df$track.name <- track.name
        track.bp.df <- track.bp.df %>% unique()
      }
    }
    if (!is.null(tracks.include)) {
      track.bp.df <- track.bp.df %>% filter(track.name %in% tracks.include)
    }
    if (cluster.tracks == T) {
      if (track.type != "bdg")
        stop("only bdg supported for clustering")
      if (!is.null(cluster.smart.cut)) {
        ## format of cluster.smart.cut: it's basically a smart.cut.df
        # data.frame(chr, start, end, regionID, int.name)
        if (! cluster.smart.cut$chr[1] %in% c("aln", "pos.out"))
          cluster.blocks <- aba.map.2.consensus(abao = abao, cluster.smart.cut)
        else
          cluster.blocks <- cluster.smart.cut
      } else {
        cluster.blocks <- NULL
      }
      tracks <- aba.cluster.tracks(track.bp.df = track.bp.df, # na.2.zero = T,
                                   cluster.blocks = cluster.blocks)
      track.bp.df <- track.bp.df %>% filter(track.name %in% tracks)
      track.bp.df$track.name <- factor(track.bp.df$track.name, levels = tracks)
    } else if (!is.null(track.order)) {
      tracks <- track.bp.df$track.name %>% unique()
      if (region.as.sort.key) {
        sort.key <- sub("\\.\\..+$", "", tracks)
      } else {
        sort.key <- tracks
      }

      if (is.function(track.order)) {
        track.order <- track.order(sort.key, abao = abao, track.bp.df = track.bp.df)
      }

      if (track.order[1] == "aln") {
        track.order <- names(abao@aln.fa)
      }

      ord <- utilsFanc::sort.by(x = sort.key, y = track.order, return.order = T, missing.method = sort.missing.method)
      tracks <- tracks[ord]
      if (order.only) {
        return(tracks)
      }
      track.bp.df$track.name <- factor(track.bp.df$track.name, levels = tracks)
    }
    if (return.matrix) {
      # if(track.type != "bdg") {
      #   stop("to return a matrix, track.type has to be bdg")
      # }
      mat <- track.bp.df[, c("pos.out", col, "track.name"), drop = F] %>%
        reshape2::acast(formula = track.name ~ pos.out, value.var = col)
      colnames(mat) <- NULL
      mat[is.na(mat)] <- 0
      return(mat)
    }
    system(paste0("mkdir -p ", out.dir))
    track.files <- track.bp.df %>% split(., f = .$track.name) %>%
      lapply(function(track) {
        track.name <- track$track.name[1]
        track <- track[, c("pos.out", col), drop = F]
        track <- track %>% dplyr::rename(end = pos.out) %>%
          dplyr::mutate(start = end -1, chr = seq.name) %>%
          dplyr::select(chr, start, end, !!as.name(col))
        if (normalize == T) {
          if (!is.null(add.constant)) {
            if (add.per.site) {
              add.constant <- add.constant * ref.length
            }
            track$score <- (track$score/(sum(track$score) + add.constant)) * scale.to
          } else {
            track$score <- (track$score/sum(track$score)) * scale.to
          }
        }
        if (track.type %in% c("bed", "rmsk")) {
          warning("add bedtools to your PATH. Otherwise you will get an error")
          track <- utilsFanc::collapse.bed.fanc(in.bed = track, bedtools.path = bedtools)
        }

        track.file <- utilsFanc::write.zip.fanc(df = track,
                                           out.file = paste0(out.dir, "/", track.name, ".", track.type),
                                           bed.shift = F)
        return(track.file)
      }) %>% unlist()
    if (push.2.igv == T) {
      igv.hubgen(files.vec = track.files, override.xml = T, override.registry = T)
    }
    return()
  } else {
    stop("currently do not support un shrunken consensus")
    # cons.seq <- abao@meta.data$consensus[[consensus.name]]$seq
    # if (is.null(cons.seq))
    #   stop("is.null(cons.seq)")
    # track.bp.df$pos.shrink <- track.bp.df$aln
  }
  # verification: pileup_test_2021-08-08.R

}

aba.align.noRepeat <- function(aba.dir, ...) {
  stop("currently not working: fa.cons.shrink somehow requires abao")
  in.fa <- Sys.glob(paste0(aba.dir, "/*in.fa"))

  if (length(in.fa) != 1) {
    stop("length(in.fa) != 1")
  }
  root <- basename(in.fa) %>% sub("_in.fa", "", .)
  noR.fa <- fa.remove.repeat(in.fa = in.fa)
  noR.aligned.fa <- utilsFanc::insert.name.before.ext(noR.fa, "aligned", "_")
  # trash <- mafft.fanc(in.fa = noR.fa, mafft = MAFFT.DEFAULT, mafft.params = MAFFT.AUTO,
  #                     aligned.fa = noR.aligned.fa)
  shrink.fa <- paste0(aba.dir, "/cons_noR/", root, "_shrink_default.fa")
  fa.cons.shrink(fa = noR.aligned.fa, out.fa =  shrink.fa, ...)
  return()
}


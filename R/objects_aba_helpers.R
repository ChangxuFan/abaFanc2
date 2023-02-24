aba.import <- function(abao) {
  if (is.character(abao)) {
    if (dir.exists(abao)) {
      abao <- paste0(abao, "/abao.Rds")
    }
    abao <- readRDS(abao)
  }
}
aba.plot.conservation <- function(abao, p.in = NULL, regions.use.regex = NULL,
                                  cons.pid = 0.9, min.pid.plot = 0.5, cons.length = 5, dist.join = 0,
                                  return.cons = F,
                                  plot.out = NULL, height = 5, width = 20) {
  mat <- abao@aln.fa %>% do.call(cbind, .)
  if (!is.null(regions.use.regex)) {
    mat <- mat %>% .[, grepl(paste0(regions.use.regex, collapse = "|"), colnames(.)), drop = F]
    if (ncol(mat) < 1) {
      stop("ncol(mat) < 1")
    }
  }
  mat <- tolower(mat)
  if (any(is.na(mat))) {
    stop("any(is.na(mat))")
  }
  
  # chars <- names(table(mat))
  pid <- apply(mat, 1, function(x) {
    tab <- table(x)
    tab <- tab[names(tab) != "-"]
    if (length(tab) == 0)
      tab <- 0
    return(max(tab))
  })
  pid <- round(pid/ncol(mat), digits = 3)
  bPid <- pid >= cons.pid
  bMat <- lapply(1:cons.length, function(i) {
    end <- length(bPid) - cons.length + i
    x <- bPid[i:end]
  }) %>% do.call(cbind, .)
  # my "smoothing by shifting method"
  bCons.site <- rowSums(bMat) == ncol(bMat)
  cons.ranges <- IRanges(start = which(bCons.site), end = which(bCons.site) + cons.length - 1) %>% 
    reduce()
  end(cons.ranges) <- pmin(end(cons.ranges) + dist.join, length(bCons.site))
  cons.ranges <- reduce(cons.ranges)
  
  if (return.cons) {
    res <- list(pid = pid, cons.ranges = cons.ranges)
    return(res)
  }
  cons.df <- data.frame(pid = pmax(pid - min.pid.plot, 0), pos = 1:length(pid))
  cons.df$bCons <- "n"
  for (i in 1:length(cons.ranges)) {
    cons.df$bCons[start(cons.ranges)[i]:end(cons.ranges)[i]] <- "c"
  }
  
  breaks <- c(min(cons.df$pid), max(cons.df$pid))
  labels <- breaks + min.pid.plot
  
  p.cons <- ggplot(cons.df, aes(x = pos, y = pid)) + 
    geom_bar(aes(fill = bCons), stat = "identity", show.legend = F) + 
    scale_x_continuous(expand = c(0.001, 0.001)) +
    scale_y_continuous(breaks = breaks, labels = labels) +
    scale_fill_manual(values = c("grey75", "red"), breaks = c("n", "c")) +
    theme_classic() + xlab("") + ylab("")+
    theme( axis.line = element_blank(), axis.ticks = element_blank())
  if (is.null(p.in)) {
    pj <- p.cons
  } else {
    p.in <- p.in + theme(axis.text.x = element_blank()) + theme(legend.position="top")
    library(cowplot)
    pj <- cowplot::plot_grid(p.in, p.cons, align = "vh", ncol = 1, rel_heights = c(10,2), axis = 'lr')
  }
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, pj, height = height, width = width, dpi = 100)
  }
  
  invisible(pj)
}

aba.plot.aln <- function(abao, use.regionID.regex = NULL,
                         tracks = NULL, p.in = NULL, gap.only = F,
                         cons.name = NULL, only.diff = F,
                         label.df = NULL,
                         to.upper = T, 
                         x.lim = NULL, rev.y = T, # ggplot2 default is: the first level of y is at bottom.
                         smart.cut.df = NULL, plot = T, 
                         plot.out = NULL,height=5, width=20,  y.axis.text.size = 8,
                         use.tile = F,
                         atcg.size = NULL, atcg.font = "mono",
                         fontface = "bold", 
                         n.wrap = NULL,
                         sub.align.out = NULL, ...) {
  # when only.diff == T: the first line will be shown in full (regardless of if it's the consensus), 
  # for all other lines, only bases different from the consensus will be shown.
  
  # label.df: it can be something like this:
  # > head(motif.df)
  #   chr start end int.name
  # 1 aln   140 151    RUNX3
  # 2 aln   201 212    TBX21
  # chr will be treated as y axis name
  
  if (!is.null(n.wrap)) {
    width = ceiling(width/n.wrap) + 1
  }
  if (is.null(atcg.size))
    atcg.size <- 2
  if (!is.null(tracks)) {
    aln.df <- lapply(tracks, function(track) {
      # browser()
      regionID <- abao@meta.data$track.df %>% filter(track.name == track) %>% pull(regionID)
      aln.df <- abao@map[[regionID]] %>% mutate(y = track)
      return(aln.df)
    }) %>% Reduce(rbind,.)
    y.levels <- tracks
  } else {
    if (!is.null(use.regionID.regex)) {
      regionIDs <- names(abao@map)
      regionIDs <- regionIDs %>% .[grepl(paste0(use.regionID.regex, collapse = "|"), .)]
      if (length(regionIDs) < 1) {
        stop("length(regionIDs) < 1")
      }
      abao@map <- abao@map[regionIDs]
    }
    aln.df <- lapply(names(abao@map), function(regionID) {
      aln.df <- abao@map[[regionID]] %>% mutate(y = regionID)
      return(aln.df)
    }) %>% Reduce(rbind, .)
    y.levels <- names(abao@map)
    
    if (!is.null(cons.name)) {
      if (is.null(abao@meta.data$consensus[[cons.name]])) {
        stop("is.null(abao@meta.data$consensus[[cons.name]])")
      }
      cons.seq <- abao@meta.data$consensus[[cons.name]]$seq %>% strsplit("") %>% unlist()
      cons.df <- data.frame(aln = 1:length(cons.seq), seq = cons.seq, ori = NA, y = cons.name)
      aln.df <- rbind(cons.df, aln.df)
      y.levels <- c(cons.name, y.levels)
    }
    
    if (only.diff == T) {
      ref <- y.levels[1]
      aln.df <- aln.df %>% group_by(aln) %>% mutate(seq = only.diff.core(seq, y, ref = ref)) %>% 
        ungroup()
    }
    if (rev.y == T)
      y.levels <- rev(y.levels)
  }
  if (to.upper == T)
    aln.df$seq <- toupper(aln.df$seq)
  
  if (!is.null(x.lim))
    aln.df <- aln.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
  if (!is.null(smart.cut.df)) {
    cut.res <- aba.smart.cut(abao = abao, track.bp.df = aln.df, smart.cut.df = smart.cut.df, ...)
    aln.df <- cut.res$out.df
    aln.df$aln <- factor(aln.df$aln, levels = cut.res$levels)
  }
  
  if (!is.null(sub.align.out)) {
    # browser()
    aln.prep <- paste0(sub.align.out, ".pre")
    seqs <- aln.df %>% split(., f = factor(.$y, levels = y.levels)) %>% 
      lapply(function(x) return(x$seq))
    seqinr::write.fasta(sequences = seqs, names = names(seqs), file.out = aln.prep)
    trash <- mafft.fanc(in.fa = aln.prep, aligned.fa = sub.align.out)
  }
  
  if (plot == T) {
    if (is.null(p.in)) {
      if (!is.null(n.wrap)) {
        tmp <- aln.df$aln %>% as.numeric()
        aln.df$row.id <- ceiling(tmp/ceiling(max(tmp)/n.wrap))
        rm(tmp)
      }
      p <- ggplot(aln.df)
      aln.df$y <- factor(aln.df$y, levels = y.levels) %>% as.numeric()
      ## key step: convert the discrete y into continuous y
    } else {
      p <- p.in
      aln.df$aln <- as.character(aln.df$aln)
    }
    
    if (gap.only == T) {
      aln.df <- aln.df %>% filter(seq %in% c("-", "N", "n"))
      # browser()
      p <- p + geom_tile(data = aln.df, 
                         mapping = aes(x=aln, y = factor(y, levels = y.levels),
                                       fill = seq),
                         inherit.aes = F)
    } else {
      if (use.tile == F) {
        p <- p + geom_text(data = aln.df, 
                           # mapping = aes(x=aln, y = factor(y, levels = y.levels),
                           #               label = seq, color = seq),
                           mapping = aes(x=aln, y = y,
                                         label = seq, color = seq),
                           inherit.aes = F, size = atcg.size, family = atcg.font, 
                           fontface = fontface, show.legend = F) +
          scale_color_manual(values = c("deepskyblue1", "orange", "red", "purple", "black" ,"grey", "grey"),
                             breaks = c("A", "C", "G", "T", "N" ,"-", "."))
      } else {
        p <- p + geom_tile(data = aln.df, mapping = aes(x=aln, y = y, fill = seq), 
                           inherit.aes = F, show.legend = F) + 
          geom_text(data = aln.df, 
                    # mapping = aes(x=aln, y = factor(y, levels = y.levels),
                    #               label = seq, color = seq),
                    mapping = aes(x=aln, y = y,
                                  label = seq),
                    inherit.aes = F, size = atcg.size, family = atcg.font, 
                    fontface = fontface, show.legend = F)
        if (is.null(p.in)) {
          p <- p + 
            scale_fill_manual(values = c("green", "orange", "red", "deepskyblue1", "grey75" ,
                                         "grey75", "grey75"),
                              breaks = c("A", "C", "G", "T", "N" ,"-", "."))
        }
      }
    if (is.null(p.in)) {
      p <- p + scale_x_continuous(expand = c(0.002, 0.002)) + 
        scale_y_continuous(breaks = 1:length(y.levels), labels = y.levels, 
                           expand = c(0.01, 0.01)) +
        xlab("") + ylab("") 
    }
    }
    if (!is.numeric(aln.df$aln) && is.null(p.in)) {
      vlines <- cut.res$cut.df.buffer$start
      p <- p + scale_x_discrete(drop = F, breaks = factor(cut.res$breaks, # %>% c(vlines), 
                                                          levels = cut.res$levels)) +
        geom_vline(xintercept = factor(vlines, levels = cut.res$levels))
      
    }
    # if (!is.null(atcg.size)) {
    #   p <- p + theme(text = element_text(size = atcg.size))
    # }
    if (is.null(p.in)) {
      p <- p + theme_classic() +
        theme( axis.line = element_blank(), axis.ticks = element_blank(), 
               axis.text.y = element_text(size = y.axis.text.size, face = "italic"))
    } else {
      p <- p + theme(axis.text.y = element_text(size = y.axis.text.size, face = "italic"))
    }
    
    if (!is.null(n.wrap)) {
      p <- p + facet_wrap(~row.id, scales = "free_x", ncol = 1) +
        theme(strip.text = element_blank())
    }
    
    if (!is.null(plot.out)) {
      dir.create(dirname(plot.out), showWarnings = F, recursive = T)
      ggsave(plot.out, p, width = width, height = height, units = "in", dpi = 100)
    }
    
    return(p)
  }
  
  return()
}

aba.plot.add.anno <- function(p, anno.df.genome = NULL, abao, bSeperate = F,
                              anno.df, anno.col, fill.col = NULL, fill.midpoint,
                              y.label = 1, rev.y = T,
                              # masking options:
                              ranges.only = NULL, ranges.flip = F, # IRanges object.
                              ##only add anno to these ranges on the raw alignment coordinates
                              ##cons related functions not yet added.
                              # plot output options:
                              plot.out = NULL, width = 20, height = 5, 
                              # rect options:
                              use.rect = F, rect.color = "black", 
                              nudge.x.left = -0.5, nudge.x.right = 0.5, 
                              nudge.y.up = 0.5, nudge.y.down = -0.5,
                              # underscore options:
                              simple.shade = F, # when this is true, simply shade the plot
                              nudge.unit = 0.3, label.nudge = 2.5,
                              alpha = 0.2) {
  # this is the master for adding motif annotation to an existing alignment plot (generated
  # by aba.plot.aln)
  # anno.df: an example would be:
  # > motif.df
  #    chr start end width strand arche.name
  # 1  aln     3  22    20      *        KLF
  # 2  aln     4  23    20      *        ETS
  # 3  aln    39  58    20      *        IRF
  ####
  # currently strand info is not used.
  # chr could be regionIDs or cons.name (more precisely the labels of y axis of the plot). Or, 
  # simply use "aln", in which case specify "y.label =", to indicate where you would like to put 
  # annotation. y = 1 will lead to annotation being placed at the bottom. It becomes the top
  # if you specify rev.y = T
  
  # for anno.df, start and end are on the coordinats of aln. 
  
  # nudges: the function originally written to add boxes around nucleotide. 
  # later changed to underscore motifs, so that it gets easier when motifs overlap with each other
  if (!is.null(anno.df.genome)) {
    # making anno.df from anno.df.genome on the fly
    if (is.character(anno.df.genome)) {
      anno.df.genome <- read.table(anno.df.genome, header = T)
    }
    anno.df <- aba.map.2.consensus(abao = abao, df = anno.df.genome)
    anno.df$chr <- anno.df$regionID
    if (!is.null(ranges.only)) {
      anno.ir <- IRanges(start = anno.df$start, end = anno.df$end)
      o <- findOverlaps(query = anno.ir, subject = ranges.only, type = "within")
      if (ranges.flip) {
        stop("ranges.flip has not been developed. Sorry...")
      }
      anno.df <- anno.df[unique(sort(queryHits(o))), ]
    }
  }
  split.id <- rep(1, nrow(anno.df))
  if (bSeperate) {
    split.id <- anno.df[, anno.col] %>% factor(.)
  }
  
  pl <- anno.df %>% split(f = split.id) %>% 
    lapply(function(anno.df) {
      g <- ggplot_build(p)
      y.levels <- g$layout$panel_params[[1]]$y$get_labels()
      ymax <- length(y.levels)
      xmax <- g$layout$panel_params[[1]]$x.sec$limits[2]
      if (simple.shade == T) {
        anno.df$chr <- anno.df$chr %>% 
          sapply(function(x) {
            id <- which(y.levels == x)
            if (length(id) > 0)
              return(id[1])
            else {
              warning(paste0(x, "not found in plot"))
              return(NA)
            }
          }) %>% as.numeric()
        anno.df <- anno.df %>% na.omit()
        anno.df <- anno.df %>% mutate(x = start + nudge.x.left, x.end = end + nudge.x.right,
                                      y = chr + nudge.y.down, y.end = chr + nudge.y.up)
        if (is.null(fill.col)) {
          fill.col <- anno.col
        } else {
          if (is.numeric(anno.df[, fill.col])) {
            anno.df[, fill.col] <- round(anno.df[, fill.col], digits = 3)
          }
          anno.df <- anno.df[(order(anno.df[, fill.col])),]
        }
        p <- p + geom_rect(data = anno.df, mapping = aes_string(xmin = "x", ymin = "y",
                                                                xmax = "x.end", ymax = "y.end",
                                                                fill = fill.col),
                           inherit.aes = T, stat = "identity", alpha = alpha)
        if (fill.col != anno.col) {
          # p <- p + scale_fill_gradient(low = "#3E3D9A", high = "white")
          g.colors <- c("white", "#E8E5F2", "#3E3D9A")
          g.values <- c(0, fill.midpoint, 1)
          n.values <- length(unique(anno.df[, fill.col]))
          if (n.values < 2) {
            g.colors[1] <- g.colors[2] <- g.colors[3]
          } 
          p <- p + scale_fill_gradientn(colors = g.colors,
                                        values = g.values)
        }
        
        # anno.df[, anno.col] <- "black"
        # p <- p + geom_segment(data = anno.df, mapping = aes_string(x = "x", y = "y", 
        #                                                            xend = "x.end", yend = "y.end",
        #                                                            color = anno.col), 
        #                       inherit.aes = F) 
        
        # anno.df.2 <- anno.df %>% mutate(start = end)
        # anno.df <- rbind(anno.df, anno.df.2) %>% arrange(chr, start)
        # p <- p  + geom_tile(data = anno.df, 
        #                    mapping = aes_string(x = "start", y = "chr", fill = anno.col),
        #                    inherit.aes = F, alpha = 0.7)
        # ggplot() + geom_tile(data = data.frame(x = c(1, 8), y = c(2,2)), mapping = aes(x = x, y = y, fill = "red"), )
        
      } else {
        trash.i <- c()
        nudge.register <- rep(0, xmax)
        for (i in 1:nrow(anno.df)) {
          df <- anno.df[i, ]
          if (df$chr != "aln") {
            y <- which(df$chr == y.levels)
            if (length(y) != 1) {
              warning(paste0(df$chr, " not found in the plot. Skipping"))
              trash.i <- c(trash.i, i)
              next
            }
          } else {
            if (rev.y == T) {
              y <- ymax - y.label + 1
            }
          }
          
          if (use.rect == T) {
            rect.df <- data.frame(x = c(df$start + nudge.x.left, df$end + nudge.x.right,
                                        df$end + nudge.x.right, df$start + nudge.x.left),
                                  y = c(y + nudge.y.down, y + nudge.y.down,
                                        y + nudge.y.up, y + nudge.y.up))
            p <- p + geom_polygon(data = rect.df, mapping = aes(x = x, y = y), 
                                  fill = NA, color = rect.color)
            anno.df$chr[i] <- y
          } else {
            # if (df$start == 207)
            #   browser()
            x.start <- df$start + nudge.x.left
            x.end <- df$end + nudge.x.right
            n.nudges <- max(nudge.register[floor(x.start):ceiling(x.end)])
            
            seg.df <- data.frame(x = x.start, x.end = x.end,
                                 y = y + 0.7 + nudge.unit * n.nudges, y.end = y + 0.7 + nudge.unit * n.nudges)
            nudge.register[floor(x.start-1):ceiling(x.end+1)] <- n.nudges + 1
            p <- p + geom_segment(data = seg.df, mapping = aes(x = x, y = y, xend = x.end, yend = y.end), 
                                  inherit.aes = F)
            anno.df$chr[i] <- seg.df$y - 1
          } 
          
          
        }
        if (length(trash.i) > 1) {
          anno.df <- anno.df[ - trash.i, ]
        }
        anno.df$chr <- as.numeric(anno.df$chr)
        anno.df <- anno.df %>% mutate(center = (start + end)/2 , chr = chr + label.nudge
        )
        p <- p + ggrepel::geom_text_repel(data = anno.df, inherit.aes = F, 
                                          mapping = aes_string(x = "center", y = "chr", label = anno.col),
                                          segment.size = 0, direction = "x", segment.alpha = 1)
      }
      p <- p + ggtitle(anno.df[, anno.col][1])
    })
  p <- scFanc::wrap.plots.fanc(plot.list = pl, n.col = 1, sub.height = height,
                               sub.width = width, plot.out = plot.out)
  if (length(pl) == 1) {
    invisible(pl[[1]])
  } else {
    invisible(p)
  }
  
}

only.diff.core <- function(seq, y, ref, symbol = ".") {
  if (sum(y == ref) != 1)
    stop("sum(y == ref) != 1")
  i <- which(y == ref)
  ref.seq <- seq[i]
  seq[seq == ref.seq] <- symbol
  seq[i] <- ref.seq
  return(seq)
}

import.bed.fanc <- function(bed, col.names = c("chr", "start", "end", "forth", "fifth", "strand"),
                            return.gr = F, no.shift = F) {
  if (grepl(".gz$", bed)) {
    bed <- sub(".gz", "",bed)
  }
  df <- read.table(bed, as.is = T, sep = "\t", quote = "", header = F)
  colnames(df) <- col.names[1:ncol(df)]
  if (no.shift != F) {
    df$start <- df$start + 1
  }
  if (return.gr == T) {
    df <- GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = T)
  }
  return(df)
}

aba.smart.cut <- function(abao, bp.vec = NULL, track.bp.df, location.col = "aln", smart.cut.df,
                          n.breaks = 10, project.breaks.to = NULL) {
  # track.bp.df is no longer required. You can use any df, as long as you supply a column
  # that encodes location
  if (is.character(smart.cut.df))
    smart.cut.df <- read.table(smart.cut.df, as.is = T, sep = "\t", quote = "", header = T)
  cut.df <- aba.map.2.consensus(abao = abao, df = smart.cut.df)
  cut.df.buffer <- cut.df %>% mutate(start = start - buffer.left, end = end + buffer.right) %>% 
    select(start, end)
  if (!is.null(bp.vec)) {
    out.df <- cut.df.buffer %>% split(f= 1:nrow(cut.df.buffer)) %>% 
      lapply(function(x) {
        out.sub <- bp.vec %>% .[ . > x$start & . < x$end]
        return(out.sub)
      }) %>% Reduce(c, .)
  } else {
    out.df <- cut.df.buffer %>% split(f= 1:nrow(cut.df.buffer)) %>% 
      lapply(function(x) {
        out.sub <- track.bp.df %>% filter(!!as.name(location.col) > x$start, !!as.name(location.col) < x$end)
        return(out.sub)
      }) %>% Reduce(rbind, .)
  }
  # browser()
  if (!is.null(n.breaks)) {
    levels <- lapply(1:nrow(cut.df.buffer), function(i) {
      return(cut.df.buffer[i, "start"]:cut.df.buffer[i, "end"])
    }) %>% Reduce(c, .) %>% sort()
    n.levels <- length(levels)
    break.interval <- ceiling(n.levels/n.breaks)
    break.pos <- which(!duplicated(floor(1:n.levels/break.interval)))
    breaks <- levels[break.pos]
    if (!is.null(project.breaks.to)) {
      breaks.proj <- aba.map.core(abao, project.breaks.to, aln = breaks)
      breaks.proj <- breaks.proj + abao@ori.df %>% dplyr::filter(regionID == project.breaks.to) %>% 
        dplyr::pull(start) - 1
    } else {
      breaks.proj <- breaks
    }
    res <- list(out.df = out.df, cut.df = cut.df, cut.df.buffer = cut.df.buffer,
                levels = levels, breaks = breaks, breaks.proj = breaks.proj, 
                project.breaks.to = project.breaks.to)
    return(res)
  } else {
    return(list(out.df = out.df, cut.df = cut.df, cut.df.buffer = cut.df.buffer))
  }
  
}

aba.map.2.consensus <- function(abao, df) {
  # expected format: chr, start, end, followed by other columns
  # regionID is the only absolutely required meta column
  if (is.character(df))
    df <- read.table(df, as.is = T, sep = "\t", quote = "", header = T)
  if ("GRanges" %in% class(df)) {
    df <- df %>% `names<-`(NULL) %>% as.data.frame()
  }
  df$regionID <- gsub("[^A-Za-z0-9]", ".", df$regionID)
  df <- df %>% split(f = 1:nrow(df)) %>% 
    lapply(function(x) {
      start.regionID <- abao@ori.df %>% filter(regionID == x$regionID) %>% pull(start)
      start.on.region <- x$start - start.regionID + 1
      end.on.region <- x$end - start.regionID + 1
      x$start <- aba.map.core(abao = abao, regionID = x$regionID, pos = start.on.region)
      x$end <- aba.map.core(abao = abao, regionID = x$regionID, pos = end.on.region)
      return(x)
    })  %>% Reduce(rbind, .)
  df$chr <- "aln"
  return(df)
}



aba.get.regionID <- function(abao, track) {
  regionID <- abao$track.df %>% filter(track.name == track) %>% pull(regionID)
  return(regionID)
}

# aba.sub.aln <- function(abao, regionID = NULL, smart.cut.df) {
#   print("miao")
#   print("miao")
# }

aba.shrink <- function(abao) {
  
}

fasta2phylip <- function(in.fasta, out.phylip = NULL) {
  if (is.null(out.phylip))
    out.phylip <- paste0(tools::file_path_sans_ext(in.fasta), ".phy")
  dat <- phylotools::read.fasta(in.fasta, clean_name = F)
  phylotools::dat2phylip(dat = dat, outfile = out.phylip)
  return(out.phylip)
}

klra.2.ly49 <- function(genes, lookup.table = "~/genomes/ly49/Klra.Ly49.Lookup.tsv") {
  if (is.character(lookup.table))
    lookup.table <- read.table(lookup.table, sep = "\t", quote = "", as.is = T, header = T)
  genes.out <- c()
  for (gene in genes) {
    if (grepl("49", gene)) {
      out <- lookup.table[tolower(lookup.table$Ly49) == tolower(gene),"Klra"]
    } else {
      out <- lookup.table[tolower(lookup.table$Klra) == tolower(gene),"Ly49"]
    }
    if (length(out) != 1)
      stop(paste0("gene ", gene, "cannot be converted"))
    genes.out <- genes.out %>% c(out)
  }
  
  return(genes.out)
}

get.exon.anderson <- function(gbk, query.df, out.fa) {
  # query.df: gene; exon. 
  gbk <- readChar(gbk, file.size(gbk))
  seq <- gbk %>% sub(".+ORIGIN", "", .) %>% sub("//.+", "", .) %>% gsub("[^atcgATCG]", "", .) %>% 
    strsplit(split = "") %>% unlist()
  gbk <- genbankr::parseGenBank(text = gbk, ret.seq = T)
  query.df$search <- paste0(query.df$gene," ", query.df$exon)
  seqs <- list()
  for (feature in gbk$FEATURES) {
    hit <- which(query.df$search == feature$note)
    if (length(hit) > 2) {
      stop("multiple hits found")
    }
    if (length(hit) > 0) {
      gene <- query.df[hit, "gene"]
      # browser()
      seqs[[gene]] <- seq[feature$start:feature$end]
    }
  }
  dir.create(dirname(out.fa), showWarnings = F, recursive = T)
  seqinr::write.fasta(sequences = seqs, names = names(seqs), file.out = out.fa)
  return(out.fa)
}


aba.track.fill.gap <- function(abao,  fill.col = NULL, track.bp.df) {
  track.bp.df <- track.bp.df %>% split(., f= factor(.$track.name, levels = unique(.$track.name))) %>% 
    lapply(function(df) {
      # regionID <- aba.get.regionID(abao, df$track.name[1])
      feature.bound <- df %>% group_by(forth) %>% summarise(start = min(aln), end = max(aln)) %>% 
        ungroup() %>% as.data.frame()
      j <- left_join(abao@map[[1]][, "aln", drop = F], df)
      # if ((j$forth %>% .[!is.na(.)] %>% .[1]) == "RMER5")
      #   browser()
      filled <- feature.bound %>% split(., f=1:nrow(.)) %>% 
        lapply(function(feature) {
          filled.sub <- j %>% filter(aln >= feature$start, aln <= feature$end)
          filled.sub$forth <- feature$forth
          filled.sub$track.name <- df$track.name[1]
          if (!is.null(fill.col))
            for (i in fill.col) {
              filled.sub[, i] <- vec.fill.gap(filled.sub[, i])
            }
          return(filled.sub)
        }) %>% Reduce(rbind, .)
      return(filled)
    }) %>% Reduce(rbind,.)
  return(track.bp.df)
}

vec.fill.gap <- function(vec) {
  for (i in 1:length(vec)) {
    if (is.na(vec[i]))
      vec[i] <- vec[i-1] + 1
  }
  return(vec)
}


aba.track.fill.gap.2 <- function(abao,  fill.col = NULL, track.bp.df) {
  track.bp.df <- track.bp.df %>% split(., f= factor(.$track.name, levels = unique(.$track.name))) %>% 
    lapply(function(df) {
      # regionID <- aba.get.regionID(abao, df$track.name[1])
      feature.bound <- df %>% group_by(forth) %>% summarise(start = min(!!as.name(fill.col)), end = max(!!as.name(fill.col))) %>% 
        ungroup() %>% as.data.frame()
      # j <- left_join(abao@map[[1]][, fill.col, drop = F], df)
      # if ((j$forth %>% .[!is.na(.)] %>% .[1]) == "RMER5")
      #   browser()
      filled <- feature.bound %>% split(., f=1:nrow(.)) %>% 
        lapply(function(feature) {
          # filled.sub <- j %>% filter(!!as.name(fill.col) >= feature$start, !!as.name(fill.col) <= feature$end)
          frame.df <- data.frame(feature$start:feature$end) %>% `colnames<-`(fill.col)
          filled.sub <- left_join(frame.df, df)
          filled.sub$forth <- feature$forth
          filled.sub$track.name <- df$track.name[1]
          return(filled.sub)
        }) %>% Reduce(rbind, .)
      return(filled)
    }) %>% Reduce(rbind,.)
  return(track.bp.df)
}

aba.cluster.tracks <- function(track.bp.df, cluster.blocks = NULL) {
  # cluster.blocks: expect the format returned by aba.map.2.consensus(): aln, start, end, int.name
  if (!is.null(cluster.blocks)) {
    if (any(!cluster.blocks$chr %in% c("aln", "pos.out"))) {
      stop('any(!cluster.blocks$chr %in% c("aln", "pos.out"))')
    }
    if (length(unique(cluster.blocks$chr)) != 1) {
      stop("length(unique(cluster.blocks$chr)) != 1")
    }
    if( !"int.name" %in% colnames(cluster.blocks)) {
      stop('!"int.name" %in% colnames(cluster.blocks)')
    }
    
    block.bp <- makeGRangesFromDataFrame(cluster.blocks, keep.extra.columns = T) %>% 
      utilsFanc::gr.bp() %>% as.data.frame() %>% dplyr::select(int.name, start)
    
    colnames(block.bp)[2] <- cluster.blocks$chr[1]
    
    signal.df <- inner_join(block.bp, track.bp.df) %>% 
      dplyr::group_by(int.name, track.name) %>% 
      dplyr::summarise(score = mean(score)) %>% 
      as.data.frame()
    
  } else {
    signal.df <- track.bp.df %>% 
      dplyr::select(pos.out, track.name, score) %>% 
      dplyr::rename(int.name = pos.out)
  }
  
  mat <- reshape2::acast(signal.df, formula = track.name ~ int.name, value.var = "score")
  mat[is.na(mat)] <- 0
  # instead of applying scaling, we directly apply pearson correlation:
  dist.o <- ClassDiscovery::distanceMatrix(t(mat), metric = "pearson")
  dend <- hclust(dist.o)
  order <- dend$labels[dend$order]
  ##### old code
  # non.all.zero <- track.bp.df %>% na.omit() %>% group_by(aln) %>% 
  #   summarise(nz = sum(score != 0) > 0) %>% filter(nz == T) %>% pull(aln)
  # 
  # track.scale.center <- track.bp.df %>% filter(aln %in% non.all.zero) %>% 
  #   group_by(track.name) %>% mutate(score = scale(score, center = T)) %>% 
  #   ungroup() %>% as.data.frame()
  # mat <- reshape2::acast(track.scale.center, formula = aln ~ track.name, value.var = "score")
  # if (na.2.zero == T) {
  #   mat[is.na(mat)] <- 0
  # }
  # tracks <- get.dendro.order(as.data.frame(mat), axis = 2)
  ### END old code
  return(order)
}

smart.cut.df.2.bed <- function(smart.cut.df, out.bed = NULL, out.track.df = NULL) {
  if (is.null(out.bed))
    out.bed <- paste0(smart.cut.df, ".bed")
  if (is.null(out.track.df))
    out.track.df <- paste0(smart.cut.df, ".track.tsv")
  smart.cut.df <- read.table(smart.cut.df, as.is = T, sep = "\t", quote = "", header = T)
  bed <- smart.cut.df %>% mutate(start = start - buffer.left, end = end + buffer.right) %>% 
    dplyr::select(chr, start ,end, int.name)
  utilsFanc::write.zip.fanc(df = bed, bed.shift = F, out.file = out.bed )
  track.df <- data.frame(regionID = smart.cut.df$regionID %>% unique(), track.name = "cut", 
                         track.file = out.bed %>% normalizePath())
  utilsFanc::write.zip.fanc(df = track.df, out.file = out.track.df, zip = F, 
                            col.names = T, row.names = F)
  return()
}

aba.loci.2.smart.cut <- function(loci, int.name, regionID, buffer.left = 0, buffer.right = 0) {
  df <- loci %>% utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = T)
  df$int.name <- int.name
  df$regionID <- regionID
  df$buffer.left <- buffer.left
  df$buffer.right <- buffer.right
  return(df)
}

aba.track.df.gen <- function(abao, regionID.include = NULL, regionID.exclude = NULL, 
                             track.file, track.name,
                             out.file = NULL) {
  # track should be just one file such as an ATAC-seq file
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  regions <- abao@ori.df$regionID
  if (!is.null(regionID.include)) {
    regions <- regions %>% .[grepl(paste0(regionID.include, collapse =  "|"), .)]
  }
  if (!is.null(regionID.exclude)) {
    regions <- regions %>% .[!grepl(paste0(regionID.exclude,collapse =  "|"), .)]
  }
  df <- data.frame(regionID = regions, track.name = track.name, track.file = track.file)
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write.table(df, out.file, sep = "\t", quote = F, col.names = T, row.names = F)
  }
  return(df)
}

aba.order.by.cat.cage <- function(track.names = NULL, abao = NULL, track.bp.df,
                                  cluster.blocks, rule.set = "r1") {
  # track.names, abao and track.bp.df has to be there as arguments, because aba.write.track expects so.
  # in fact we only need track.bp.df to figure out the order.
  # blocks: expect the format returned by aba.map.2.consensus(): aln, start, end, int.name
  
  ## the same as aba.cluster.tracks:
  if (any(!cluster.blocks$chr %in% c("aln", "pos.out"))) {
    stop('any(!cluster.blocks$chr %in% c("aln", "pos.out"))')
  }
  if (length(unique(cluster.blocks$chr)) != 1) {
    stop("length(unique(cluster.blocks$chr)) != 1")
  }
  if( !"int.name" %in% colnames(cluster.blocks)) {
    stop('!"int.name" %in% colnames(cluster.blocks)')
  }
  
  block.bp <- makeGRangesFromDataFrame(cluster.blocks, keep.extra.columns = T) %>% 
    utilsFanc::gr.bp() %>% as.data.frame() %>% dplyr::select(int.name, start)
  
  colnames(block.bp)[2] <- cluster.blocks$chr[1]
  
  if (rule.set == "r1") {
    # r1 functions by max
    signal.df <- inner_join(block.bp, track.bp.df) %>% 
      dplyr::group_by(int.name, track.name) %>% 
      dplyr::summarise(score = max(abs(score))) %>% 
      as.data.frame() %>% 
      reshape2::dcast(track.name ~ int.name, value.var = "score")
    if (any(!c("exon2", "exon1", "RMER5") %in% colnames(signal.df))) {
      stop('any(!c("exon2", "exon1", "RMER5") %in% colnames(signal.df))')
    }    
    cat.df <- signal.df %>% split(., 1:nrow(.)) %>% 
      lapply(function(df) {
        track.name <- df$track.name
        df$track.name <- NULL
        df[is.na(df)] <- 0
        cat.df <- data.frame(track.name = track.name, cat = NA)
        if (df$RMER5 > 0.3 * max(df)) {
          cat.df$cat <- 1
        } else if (df$exon2 > 0.3 * max(df)) {
          cat.df$cat <- 3
        } else {
          cat.df$cat <- 2
        }
        return(cat.df)
      }) %>% Reduce(rbind, .) %>% 
      dplyr::arrange(cat, track.name)
  } else {
    stop("only rule set r1 has been developed")
  }
  res <- cat.df$track.name
  if (!is.null(track.names)) {
    if(!grepl("\\.\\.", track.names[1])) {
      res <- res %>% sub("\\.\\..+$", "", .)
    }
  }
  return(res)
}

aba.shift <- function(abao, shift.start, shift.end, 
                      buffer.start = 0, buffer.end = 0, ref.region,
                      regionIDs.include = NULL, regionIDs.exclude = NULL,
                      new.dir, threads = 1, ...) {
  aba.df <- aba.subset(abao = abao, regionIDs = regionIDs.include, 
                       regionIDs.exclude = regionIDs.exclude, ori.df.only = T)
  aba.df$start <- pmax(aba.df$start + shift.start + buffer.start, 0)
  aba.df$end <- pmax(aba.df$end + shift.end + buffer.end, aba.df$start)
  zero <- which(aba.df$end - aba.df$start == 0)
  if (length(zero) > 0) {
    stop("the following regions had no width after shifting: \n" %>% 
           paste0(paste0(aba.df$regionID[zero]), collapse = "; "))
  }
  if (abs(buffer.start) + abs(buffer.end) > 0) {
    if (buffer.start > 0) {
      stop("buffer.start must <=0 to make the region larger")
    }
    if (buffer.end < 0) {
      stop("buffer.end must >=0 to make the region larger")
    }
    smart.cut <- aba.df %>% dplyr::filter(regionID == ref.region) %>% 
      dplyr::mutate(start = start - buffer.start, end = end - buffer.end)
    
    tmp.dir <- paste0(sub("/+$", "", new.dir), "_buffer/")
    abao.tmp <- aba.create(aba.df = aba.df, df.is.bed = F, threads = threads, 
                           work.dir = tmp.dir, ...)

    abao <- aba.subset(abao.tmp, smart.cut.df = smart.cut, threads = threads, 
                       new.work.dir = new.dir)
  } else {
    abao <- aba.create(aba.df = aba.df, df.is.bed = F, threads = threads, 
                       work.dir = new.dir, ...)
  }
  invisible(abao)
}
mb.tree.shorten.branch <- function(tree.in, tree.out = NULL, prob.percent.cutoff = 85) {
  if (is.null(tree.out)) {
    tree.out <- tools::file_path_sans_ext(tree.in) %>% 
      paste0("_shorten", prob.percent.cutoff, ".tre")
  }
  tree <- readLines(tree.in)
  bKey <- grepl("tree.+con", tree)
  if (sum(bKey) != 1) {
    stop("sum(bKey) != 1")
  }
  key <- tree[bKey]
  keys <- key %>% strsplit(":") %>% unlist()
  len <- length(keys) - 1
  
  prob.pcts <- sub("^.+prob.percent.", "", keys[1:len]) %>% 
    sub("prob.+$", "", .) %>% stringr::str_extract("\\d+") %>% as.numeric()
  
  bShorten <- prob.pcts < prob.percent.cutoff
  bShorten <- c(F, bShorten)
  
  keys[bShorten] <- sub("e.\\d+", "e-10",keys[bShorten])
  
  key <- paste0(keys, collapse = ":")
  tree[bKey] <- key
  
  dir.create(dirname(tree.out), showWarnings = F, recursive = T)
  write(tree, tree.out, sep = "\n")
  return(tree.out)
}

treedata.get.label <- function(td) {
  tb <- as_tibble(td)
  return(tb$label %>% .[!is.na(.)])
}

ggtree.fanc <- function(tree, read.fun = treeio::read.tree, out.file = NULL, 
                        layout = "daylight", bUse.tiplab = F,
                        group.map = NULL, color.map = NULL, tip.lab.size = 4, 
                        tree.width = 7, tree.height = 5) {
  if (is.character(tree)) {
    tree <- read.fun(tree)
  }
  tree <- treeio::as.treedata(tree)
  if (!is.null(group.map)) {
    if (is.character(group.map)) {
      group.map <- read.table(group.map, header = T, comment.char = "", quote = "", sep = "\t") %>% 
        .[, c("label", "group")]
    }
    group.map <- group.map %>% dplyr::filter(label %in% treedata.get.label(tree))
    if (nrow(group.map) < 1) {
      stop("nrow(group.map) < 1")
    }
    tree <- full_join(tree, group.map, by = "label")
  } else {
    tmp <- tree %>% as_tibble()
    tmp.map <- data.frame(label = tmp$label %>% .[!is.na(.)] %>% .[.!= ""] %>% 
                            .[is.na(suppressWarnings(as.numeric(.)))], 
                          group = "ungrouped")
    tree <- full_join(tree, tmp.map, by = "label")
  }
  
  if (!is.null(color.map)) {
    if (is.character(color.map)) {
      color.map <- read.table(color.map, header = T, comment.char = "", quote = "", sep = "\t") %>% 
        .[, c("group", "color")]
    }
    colors <- color.map$color
    names(colors) <- color.map$group
  } else {
    colors <- c("ungrouped" = "black")
  }

  # if (layout %in% c("daylight", "equal_angle")) {
  #   t <- tree %>% as_tibble() 
  #   t$label[sapply(t$label, is.numeric)] <- NA
  #   ggtree(as.treedata(t), layout = layout)
  # }
  
  p <- ggtree(tree, layout = layout)
  
  if (layout == "rectangular" || bUse.tiplab ) {
    p <- p + geom_tiplab(aes(color = group), size = tip.lab.size)
  } else {
    p <- p + ggrepel::geom_text_repel(aes(label = label, color = group), size = tip.lab.size)
  }
  p <- p + scale_color_manual(values = colors) + coord_cartesian(clip="off")
  
  y.lims <- layer_scales(p)$y$range$range %>% sort()
  x.lims <- layer_scales(p)$x$range$range %>% sort()
  p <- p + geom_treescale(x = x.lims[1], y = y.lims[2], offset = 0.05 * (y.lims[2] - y.lims[1]))
  p <- p + scale_x_continuous(expand = expansion(mult = 0.25)) +
    theme(aspect.ratio = 1)
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    ggsave(out.file, p, width = tree.width, height = tree.height, dpi = 100, unit = "in")
  }
  return(p)
}

tree.from.mat <- function(mat, dist.method = "euclidean", per.row = F,
                          signal.outlier.quantile = 0.02, force = F,
                          norm.factors = NULL,
                          bLog2p1 = T, bLog.before.norm = F, log.addon = 1,
                          bLog.dist.mat = F,
                          tree.method = "upgma",
                          scale.to = 100, bStandardize = F,
                          bSmooth = F, span = 0.1) {
  
  # mat: each row is a Ly49, each column is a feature (such as the ATAC-seq signal at a nucleotide)
  mat.ori <- mat
  if (! is.null(signal.outlier.quantile)) {
    if (ncol(mat) < 2) {
      warning(paste0("each Ly49 is represented by 1 number. filtering by quantile can be a bad idea\n",
                     "therefore not performing outlier removal. use force = T to override."))
    } else {
      warning("removing signal outliers!!")
      remove.outlier <- function(x) {
        low.bound <- quantile(x, signal.outlier.quantile)
        high.bound <- quantile(x, 1 - signal.outlier.quantile)
        x[x < low.bound] <- low.bound
        x[x > high.bound] <- high.bound
        return(x)
      }
      
      if (per.row) {
        mat <- apply(mat, 1, remove.outlier) %>% t()
      } else {
        mat <- remove.outlier(mat)
      }
    }
  } else {
    print("not removing signal outliers!!")
  }
  if (bLog2p1 && bLog.before.norm) {
    print("transforming before norm")
    mat <- log2(mat + log.addon)
  }
  if (!is.null(scale.to) && !is.na(scale.to)) {
    # scale.f <- function(x) {
    #   x <- round(x * scale.to * length(x) / sum(x), 3)
    # }
    # 
    # if (per.row) {
    #   mat <- apply(mat, 1, scale.f) %>% t()
    # } else {
    #   mat <- scale.f(mat)
    # }
    
    ## upgrade: double scaling: 
    # first round scaling: scale rowsums to the same mean across datasets.
    # bk <- mat
    mat[rowSums(mat) == 0, 1] <- 1
    if (per.row) {
      rsums <- rep(scale.to, nrow(mat))
    } else {
      rsums <- rowSums(mat) * scale.to / mean(rowSums(mat))
    }
    # second round scaling: scale each row to the rowsums you got from the first round (rsums)
    rnames <- rownames(mat)
    mat <- diag(rsums) %*% diag(1/rowSums(mat)) %*% mat
    rownames(mat) <- rnames
    # rationale for double scaling: within each row, we want to scale to the sum, because we want 
    # a relatively steady matrix size when you have either fewer or more columns
    # between rows, we scale rowsums to the same mean across datasets, so that the scaling is not 
    # affected by adding or removing strains.
  } else {
    warning("no scaling used!")
  }
  
  if (bLog2p1 && !bLog.before.norm) {
    print("transforming after norm")
    mat <- log2(mat + log.addon)
  }
  
  if (bStandardize) {
    print("standardizing")
    mat <- scale(mat)
  }
  
  if (bSmooth) {
    mat_unsmoothed <- mat
    mat <- apply(mat, 1, function(signal) {
      df <- data.frame(x = 1:ncol(mat), y = signal)
      l <- loess(y ~ x, df, degree = 2, span = span)
      smoothed <- predict(l, df)
      return(smoothed)
    }) %>% t()
    colnames(mat) <- NULL
  } else {
    mat_unsmoothed <- NULL
  }
  if (dist.method %in% c("pearson", "spearman", "sqrt pearson", "absolute pearson")) {
    dist.o <- ClassDiscovery::distanceMatrix(t(mat), metric = dist.method)
  } else {
    dist.o <- dist(mat, method = dist.method)
  }
  
  if (bLog.dist.mat) {
    dist.o <- log2(dist.o + log.addon)
  }
  if (tree.method == "upgma") {
    print("using upgma tree")
    tree <- phangorn::upgma(dist.o)
  } else if (tree.method == "fastme") {
    print("using fastme tree")
    tree <- ape::fastme.ols(dist.o, nni = T)
  } else {
    stop("only upgma and fastme trees are supported!")
  }
  
  res <- list(mat = mat, dist.o = dist.o, tree = tree,
              mat_unsmoothed = mat_unsmoothed,
              tree.method = tree.method, mat.ori = mat.ori)
  return(res)
}

tree.build.viz.core <- function(mat, out.dir, root.name,
                                signal.outlier.quantile = 0.02,
                                bBy.row = F, bCounts.only = F, bMap18.mask = F,
                                bLog2p1 = T, bLog.before.norm = F, log.addon = 1,
                                bLog.dist.mat = F,
                                scale.factor = 100, bStandardize = F,
                                dist.method = "euclidean", tree.method = "upgma",
                                bSmooth = F, span = 0.1,
                                group.map = NULL, colors,
                                hm.force.row.order = T, hm.cell.font.size = 10,
                                tip.lab.size = 4, tree.width = 7, tree.height = 7) {
  by.row <- ifelse(bBy.row, "_row", "_full")
  smooth <- ifelse(bSmooth, paste0("_", span), "")
  if (bCounts.only && !bBy.row) {
    if (bMap18.mask) {
      rnames <- rownames(mat)
      mat <- cbind(rowSums(mat[, 1:250]), rowSums(mat[, 251:500]))
      rownames(mat) <- rnames
      rm(rnames)
    } else {
      rnames <- rownames(mat)
      mat <- matrix(rowSums(mat), ncol = 1)
      rownames(mat) <- rnames
      rm(rnames)
    }
    
  }
  tree <- tree.from.mat(mat, per.row = bBy.row, signal.outlier.quantile = signal.outlier.quantile,
                        bLog2p1 = bLog2p1, bLog.before.norm = bLog.before.norm, log.addon = log.addon,
                        bLog.dist.mat = bLog.dist.mat,
                        scale.to = scale.factor, bStandardize = bStandardize,
                        tree.method = tree.method,
                        bSmooth = bSmooth, 
                        span = span, 
                        dist.method = dist.method)
  if (tree.method == "fastme") {
    tree$tree <- phytools::midpoint.root(tree = tree$tree)
  }
  tree$tree <- as.treedata(tree$tree)
  if (!is.null(group.map)) {
    group.map <- group.map %>% dplyr::filter(label %in% treedata.get.label(tree$tree))
    if (nrow(group.map) < 1) {
      stop("nrow(group.map) < 1")
    }
    tree$tree <- full_join(tree$tree, group.map, by = "label")
  } else {
    tmp <- tree$tree %>% as_tibble()
    tmp.map <- data.frame(label = tmp$label %>% .[!is.na(.)], group = "ungrouped")
    tree$tree <- full_join(tree$tree, tmp.map, by = "label")
  }
  
  # cell_fun <- function(j, i, x, y, width, height, fill) {
  #   grid.text(round(mat[i, j], digits = 2), x, y, 
  #             gp = gpar(fontsize = hm.cell.font.size))
  # }
  cell_fun <- NULL
  
  cell_fun2 <- function(j, i, x, y, width, height, fill) {
    grid.text(round(as.matrix(tree$dist.o)[i, j], digits = 2), x, y, 
              gp = gpar(fontsize = hm.cell.font.size))
  }
  
  if (hm.force.row.order) {
    p <- ComplexHeatmap::Heatmap(tree$mat, cluster_rows = F, cluster_columns = F, cell_fun = cell_fun,
                                 row_order = attr(tree$dist.o, "Labels"))
    p2 <- ComplexHeatmap::Heatmap(as.matrix(tree$dist.o), cluster_rows = F, cluster_columns = F,
                                  cell_fun =  cell_fun2,
                                  row_order = attr(tree$dist.o, "Labels"),
                                  column_order = attr(tree$dist.o, "Labels"))
  } else {
    p <- ComplexHeatmap::Heatmap(tree$mat, cluster_rows = T, cluster_columns = F,
                                 cell_fun = cell_fun,
                                 clustering_distance_rows = "euclidean")
    p2 <- ComplexHeatmap::Heatmap(as.matrix(tree$dist.o), cluster_rows = T, cluster_columns = T,
                                  cell_fun = cell_fun2)
  }
  hm.file <- paste0(out.dir, "/", root.name, by.row, smooth, "_mat.png")
  scFanc::save.base.plot(p = p, file = hm.file, width = 1400, height = 700)
  
  scFanc::save.base.plot(p = p2, width = 900, height = 700,
                         file = paste0(out.dir, "/", root.name, by.row, smooth, "_dist_mat.png"))
  
  
  if (bSmooth) {
    pl <- lapply(rownames(tree$mat), function(gene) {
      df <- data.frame(pos = 1:ncol(tree$mat), 
                       un_smoothed = tree$mat_unsmoothed[gene, ], smoothed = tree$mat[gene, ]) %>% 
        reshape2::melt(id.vars = "pos", value.name = "signal", variable.name = "type")
      p <- ggplot(df, aes(x = pos,  y = signal)) + 
        geom_line(aes(group = type, linetype = type, color = type)) +
        ggtitle(gene)
      return(p)
    })
    trash <- scFanc::wrap.plots.fanc(pl, plot.out = sub("_mat", "_gaze", hm.file), sub.width = 7)
  }
  lapply(c("phylogram", "unrooted"), function(tree.type) {
    tree.out <- paste0(out.dir, "/", root.name, by.row, smooth, "_", tree.type, ".png")
    dir.create(dirname(tree.out), showWarnings = F, recursive = T)
    # if (!is.null(color.map)) {
    #   tip.colors <- data.frame(tip = tree$tree$tip.label) %>% 
    #     left_join(color.map) %>% .[!duplicated(.$tip),] %>% pull(color)
    # } else {
    #   tip.colors <- "black"
    # }
    # browser()
    # 
    # sbp <- scFanc::save.base.plot
    # environment(sbp) <- environment()
    # trash <- sbp(exp = ape::plot.phylo(tree$tree, tree.type, tip.color = tip.colors), 
    #              file = tree.out)
    
    if (tree.type == "phylogram") 
      layout = "rectangular"
    else if (tree.type == "unrooted")
      layout = "daylight"
    else
      stop("tree.type not recognized")
    
    p <- ggtree(tree$tree, layout = layout)
    if (tree.type == "phylogram") {
      p <- p + geom_tiplab(aes(color = group), size = tip.lab.size)
    } else {
      p <- p + geom_text_repel(aes(label = label, color = group), size = tip.lab.size)
    }
    p <- p + scale_color_manual(values = colors) + coord_cartesian(clip="off")
    
    y.lims <- layer_scales(p)$y$range$range %>% sort()
    x.lims <- layer_scales(p)$x$range$range %>% sort()
    p <- p + geom_treescale(x = x.lims[1], y = y.lims[2])
    p <- p + scale_x_continuous(expand = expansion(mult = 0.25))
    p <- p + ggsave(tree.out, width = tree.width, height = tree.height, dpi = 100, unit = "in")
    return()
  })
  return(tree)
}

epi.tree.pipe.2 <- function(abao, signal.files, bMap18.mask = F, signal.outlier.quantile = 0.02,
                            bCounts.only = F, bLog2p1 = T, bLog.before.norm = F, log.addon = 1,
                            bLog.dist.mat = F,
                            scale.factor = 100, bStandardize = F, norm.factors = NULL,
                            dist.method = "euclidean", tree.method = "upgma",
                          bSmooth = F, span = 0.1, expand.smooth.1side = NULL,
                          bCat.samples = T, dist.cor.groups = NULL, sampleName.sep = "_",
                          group.map = NULL, color.map = NULL,
                          hm.force.row.order = T,
                          tip.lab.size = 4, bSimplify.label = T, hm.cell.font.size = 10,
                          out.dir, threads = 1) {
  # color.map format: data.frame(group = , color = )
  # group.map format: data.frame(label = "Ly49h", group = "acti_1")
  if (!is.null(expand.smooth.1side)) {
    bSmooth <- F
  }
  if (is.null(names(signal.files))) {
    stop("signal.files must be named. These names will be used as rootnames for output")
  }
  
  if (!is.null(color.map) && is.character(color.map)) {
    color.map <- read.table(color.map, header = T, comment.char = "", quote = "", sep = "\t") %>% 
      .[, c("group", "color")]
    colors <- color.map$color
    names(colors) <- color.map$group
  }
  if (is.null(color.map)) {
    colors <- c("ungrouped" = "black")
  }
  
  if (!is.null(group.map) && is.character(group.map)) {
    group.map <- read.table(group.map, header = T, comment.char = "", quote = "", sep = "\t") %>% 
      .[, c("label", "group")]
    if (bSimplify.label) {
      group.map$label <- sub("^.*Ly49", "", group.map$label)
    }
  }
  
  smooth <- ifelse(bSmooth, paste0("_", span), "")
  # abao <- aba.add.consensus(abao, shrink.all = T, force = T)
  if (!is.null(norm.factors)) {
    print("using norm.factors, scale.factor is automatically turned off")
    scale.factor <- NULL
  }
  
  forest <- utilsFanc::safelapply(names(signal.files), function(root.name)  {
    out.dir <- paste0(out.dir, "/", root.name, "/")
    dir.create(out.dir, showWarnings = F, recursive = T)
    signal.file <- signal.files[[root.name]]
    abao <- aba.add.track(abao, track.file = signal.file, track.type = "bdg", 
                          smooth.win.1side = expand.smooth.1side)
    mat <- aba.write.track(abao, consensus.name = "cons_N_0.5_A", 
                           track.name = "bdg", track.type = "bdg", return.matrix = T)
    rownames(mat) <- rownames(mat) %>% sub("\\.\\..+$", "", .)
    if (bMap18.mask) {
      mat <- mat[, c(1:250, (ncol(mat) - 249):ncol(mat))]
    }
    if (bSimplify.label) {
      rownames(mat) <- sub("^.*Ly49", "", rownames(mat))
    }
    if (!is.null(norm.factors)) {
      norm.factor <- norm.factors[root.name]
      if (is.na(norm.factor)) {
        stop("is.na(norm.factor)")
      }
      mat <- mat / norm.factor
    }
    
    trees <- lapply(c(F, T), function(bBy.row) {
      tree <- tree.build.viz.core(mat = mat, out.dir = out.dir, bBy.row = bBy.row,
                                  signal.outlier.quantile = signal.outlier.quantile,
                                  root.name = root.name, bCounts.only = bCounts.only, 
                                  bMap18.mask = bMap18.mask, bLog2p1 = bLog2p1,
                                  bLog.before.norm = bLog.before.norm, log.addon = log.addon, 
                                  bLog.dist.mat = bLog.dist.mat, scale.factor = scale.factor, 
                                  bStandardize = bStandardize, dist.method = dist.method, 
                                  tree.method = tree.method, bSmooth = bSmooth, span = span,
                                  group.map = group.map, colors = colors,
                                  hm.force.row.order = hm.force.row.order, hm.cell.font.size = hm.cell.font.size,
                                  tip.lab.size = tip.lab.size)
      
      return(tree)
    })
    names(trees) <- c("full", "row")
    saveRDS(trees, paste0(out.dir, "/", root.name, "_trees.Rds"))
    return(trees)
  }, threads = threads)
  names(forest) <- names(signal.files)
  if (bSmooth) {
    span.flag <- paste0("_", span)
  } else {
    span.flag <- ""
  }
  if (bCounts.only && !bMap18.mask) {
    scaled.mat <- lapply(forest, function(x) return(x$full$mat)) %>% Reduce(cbind, .)
    colnames(scaled.mat) <- names(forest)
    cell_fun <- function(j, i, x, y, width, height, fill) {
      grid.text(round(scaled.mat[i, j], digits = 2), x, y, 
                gp = gpar(fontsize = hm.cell.font.size))
    }
    p <- ComplexHeatmap::Heatmap(scaled.mat, cell_fun = cell_fun)
    scFanc::save.base.plot(p = p, file = paste0(out.dir, "/scaled_mat.png"))
    saveRDS(scaled.mat, paste0(out.dir, "/scaled.mat.Rds"))
  }
  ## build a tree with all samples catted into one matrix
  if (bCat.samples) {
    cat.mat <- lapply(names(forest), function(root.name) {
      x <- forest[[root.name]]
      mat <- x$full$mat.ori
      rownames(mat) <- paste0(rownames(mat), "_", root.name)
      return(mat)
    }) %>% Reduce(rbind, .)
    dist.ref <- Biostrings::readDNAMultipleAlignment(paste0(abao@work.dir, "/cons/cons_N_0.5_A/", 
                                                            abao@meta.data$fa.root.name ,"cons_N_0.5_A_shrink.fa")) %>% 
      Biostrings::as.matrix() %>% ape::as.DNAbin() %>% ape::dist.dna(as.matrix = T, model = "TN93")
    
    cat.trees <- lapply(c(T,F), function(bBy.row) {
      by.row <- ifelse(bBy.row, "row", "full")
      cat.group.map <- lapply(names(signal.files), function(root.name) {
        df <- group.map
        df$label <- paste0(df$label, "_", root.name)
        return(df)
      }) %>% Reduce(rbind, .)
      
      cat.tree <- tree.build.viz.core(mat = cat.mat, out.dir = paste0(out.dir, "/cat_samples/"),
                                      root.name = "cat_samples", 
                                      signal.outlier.quantile = signal.outlier.quantile,
                                      bBy.row = bBy.row, bCounts.only = bCounts.only,
                                      bMap18.mask = bMap18.mask, bLog2p1 = bLog2p1,
                                      bLog.before.norm = bLog.before.norm, log.addon = log.addon, 
                                      bLog.dist.mat = bLog.dist.mat, scale.factor = scale.factor, 
                                      bStandardize = bStandardize, dist.method = dist.method, 
                                      tree.method = tree.method, bSmooth = bSmooth, span = span,
                                      group.map = cat.group.map, colors = colors,
                                      hm.force.row.order = F, tip.lab.size = tip.lab.size,
                                      tree.width = 14, tree.height = 14)
      dist.cor(dist.mat.query.list = as.matrix(cat.tree$dist.o), dist.mat.ref = dist.ref, 
               sampleName.sep = sampleName.sep, ref.groups = dist.cor.groups, 
               Ly49.simplify = bSimplify.label, jmat.scale.to = 1, out.dir = paste0(out.dir, "/cat_samples/"),
               root.name = by.row)
    })
    names(cat.trees) <- c("full", "row")
    saveRDS(cat.trees, paste0(out.dir, "/cat_trees.Rds"))
  }
  
  plot.types <- outer(paste0(c("full", "row"), span.flag), c("mat", "dist_mat", "phylogram", "unrooted"),
                      paste, sep = "_") %>% as.vector()
  utilsFanc::safelapply(plot.types, function(plot.type) {
    plots <- paste0(out.dir, "/", names(signal.files), "/", 
                    names(signal.files), "_", plot.type,  ".png")
    utilsFanc::png.wrap(pngs = plots, width = 10, height = 10, add.title = T,
                        plot.out = paste0(out.dir, "/combined_plots/", plot.type, ".png"))
  }, threads = threads)
  saveRDS(forest, paste0(out.dir, "/forest.Rds"))
  return(forest)
}
get.tip.color <- function(tree, color.map) {
  if (is.character(color.map)) {
    color.map <- read.table(color.map, header = T, comment.char = "", quote = "", sep = "\t")
  }
  colors <- data.frame(tip = tree$tip.label) %>% 
    left_join(color.map) %>% .[!duplicated(.$tip),] %>% pull(color)
  return(colors)
}

dist.cor <- function(dist.mat.query.list, dist.mat.ref, sampleName.sep = NULL, ref.groups = NULL,
                     Ly49.simplify = T,
                     jmat.scale.to = 10, 
                     out.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(out.dir)
  if (!is.null(sampleName.sep)) {
    t.mat <- dist.mat.query.list
    if (!is.matrix(t.mat) || !identical(t.mat, t(t.mat))) {
      stop("!is.matrix(t.mat) || !identical(t.mat, t(t.mat))")
    }
    seqs <- colnames(t.mat)
    groups <- strsplit(seqs, sampleName.sep) %>% lapply(function(x) return(x[length(x)])) %>% 
      unlist()
    groups <- factor(groups, levels = unique(groups))
    dist.mat.query.list <- seqs %>% split(f = groups) %>% 
      lapply(function(x) {
        res <- t.mat[x, x]
        rownames(res) <- colnames(res) <- sub(paste0(sampleName.sep, ".*$"), "", rownames(res))
        return(res)
      })
    names(dist.mat.query.list) <- as.character(names(dist.mat.query.list))
  }
  if (is.null(names(dist.mat.query.list))) {
    stop("is.null(names((dist.mat.query.list))")
  }
  dist.mat.ref <- list(ref = dist.mat.ref)
  mats <- c(dist.mat.ref, dist.mat.query.list)
  rm(dist.mat.ref)
  
  mat.names <- names(mats)
  mats <- lapply(names(mats), function(mat.name) {
    mat <- mats[[mat.name]]
    if (!identical(mat, t(mat))) {
      stop(paste0(mat.name, " is not a symmetrical matrix"))
    }
    if (Ly49.simplify) {
      rownames(mat) <- sub("^.*Ly49", "", rownames(mat))
      colnames(mat) <- sub("^.*Ly49", "", colnames(mat))
    }
    return(mat)
  })
  names(mats) <- mat.names; rm(mat.names)
  
  if (is.null(ref.groups)) {
    ref.groups <- list(group1 = colnames(mats[[1]]))
  } else {
    if (!is.list(ref.groups)) {
      group.names <- basename(ref.groups)
      ref.groups <- lapply(ref.groups, readLines)
      names(ref.groups) <- group.names %>% tools::file_path_sans_ext()
    }
    ref.groups$all <- unlist(ref.groups) %>% unique()
  }
  if (Ly49.simplify) {
    ref.groups <- lapply(ref.groups, function(x) return(sub("^.*Ly49", "", x)))
  }
  pl <- lapply(names(ref.groups), function(group.name) {
    to.use <- ref.groups[[group.name]]
    
    regions <- lapply(mats, colnames) %>% Reduce(intersect, .)
    if (length(regions) < 2) {
      warning(paste0("ref.group: ",group.name, ": length(regions) < 2"))
      return()
    }
    
    to.use <- intersect(to.use, regions)
    if (length(to.use) < 2) {
      warning(paste0("ref.group: ",group.name, 
                     ": length(to.use) < 2 after intersecting regions in matrices"))
      return()
    }
    # lapply(names(mats), function(mat.name) {
    #   utilsFanc::check.intersect(to.use, paste0("genes in ref.group: ", group.name),
    #                              colnames(mats[[mat.name]]), paste0("genes in matrix: ", mat.name))
    # })
    
    mats <- lapply(mats, function(x) return(x[to.use, to.use]))
    
    # qc using heatmaps:
    hm.ref <- ComplexHeatmap::Heatmap(mats$ref, name = "ref")
    hms <- lapply(names(mats), function(mat.name) {
      mat <- mats[[mat.name]]
      hm <- ComplexHeatmap::Heatmap(mat, column_order = suppressWarnings(column_order(hm.ref)),
                                    row_order = suppressWarnings(row_order(hm.ref)), name = mat.name)
      return(hm)
    })
    hms[[1]] <- hm.ref
    hm.big <- hms[[1]]
    for (i in 2:length(hms)) {
      hm.big <- hm.big + hms[[i]]
    }
    scFanc::save.base.plot(p = hm.big, file = paste0(out.dir, "/", root.name, "_", group.name, "_mats.png"), 
                           height = 50 * ncol(mats[[1]]), width = 50 * ncol(mats[[1]]) * length(hms) + 200)
    # another qc: pileup all the entries from all samples:
    jmat <- lapply(mats, function(mat) {
      vec <- mat[upper.tri(mat)]
      return(vec)
    }) %>% Reduce(rbind, .)
    rownames(jmat) <- names(mats)
    
    name.mat <- outer(to.use, to.use, paste, sep = "_")
    nvec <- name.mat[upper.tri(name.mat)]
    colnames(jmat) <- nvec
    jmat.norm <- diag(rep(jmat.scale.to, nrow(jmat))) %*% diag(1/rowMeans(jmat)) %*% jmat
    rownames(jmat.norm) <- rownames(jmat)
    hm.jmat <- ComplexHeatmap::Heatmap(jmat.norm, cluster_rows = T,
                                       cluster_columns = T)
    
    scFanc::save.base.plot(p = hm.jmat, file = paste0(out.dir, "/", root.name, "_", group.name, "_jmat.png"),
                           height = 50 * nrow(jmat.norm), width = 30 * ncol(jmat.norm))
    
    ## correlation plot:
    df.cor <- data.frame(ref = jmat["ref", ], 
                         query.mean = colMeans(jmat[-1, ]),
                         pair = colnames(jmat))
    p.cor <- ggplot(df.cor, aes(x = ref, y = query.mean, label = pair)) + 
      geom_point() + 
      ggrepel::geom_text_repel(color = "blue") +
      ggpubr::stat_cor(method = "pearson") +
      ggtitle(group.name) + 
      theme(aspect.ratio = 1) +
      ggsave(paste0(out.dir, "/", root.name, "_", group.name, "_cor.pdf"), width = 5.5, height = 5, units = "in")
    df.bar <- reshape2::melt(jmat.norm, varnames = c("id", "pair"), value.name = "dist") %>% 
      scFanc::factor2character.fanc()
    df.bar$type <- df.bar$id
    df.bar$type[df.bar$type != "ref"] <- "query"
    df.bar.ref <- df.bar %>% dplyr::filter(type == "ref")
    df.bar.query <- df.bar %>% dplyr::filter(type == "query")
    df.bar$pair <- factor(df.bar$pair, levels = df.bar.ref$pair[order(df.bar.ref$dist)])
    df.bar.query$pair <- factor(df.bar.query$pair, levels = df.bar.ref$pair[order(df.bar.ref$dist)])
    # ggpubr::ggbarplot(df.bar, x = "pair", y = "dist", fill = "type", add = "mean_sd", 
    #                   position = position_dodge(0.9), alpha = 0.5)
    p.rank <- ggpubr::ggbarplot(df.bar.query, x = "pair", y = "dist", 
                                fill = "seagreen", alpha = 0.3, add = "mean_sd") + 
      geom_point(data = df.bar.ref) +
      ggtitle(group.name) +
      scale_x_discrete(guide = guide_axis(angle = 90)) + 
      ggsave(paste0(out.dir, "/", root.name, "_", group.name, "_rank.pdf"),
             width = 0.5 * length(df.bar.ref$pair %>% unique()), height = 5, units = "in",
             limitsize = F)
    res <- list(p.cor = p.cor, p.rank = p.rank)
    return(res)
  })
  pl <- pl[!is.null(pl)]
  if (length(pl) < 1) {
    stop("length(pl) < 1")
  }
  names(pl) <- names(ref.groups)
  pl.cor <- lapply(pl, function(x) return(x$p.cor))
  trash <- scFanc::wrap.plots.fanc(pl.cor, plot.out = paste0(out.dir, "/", root.name, "_allGroups_cor.png"), 
                          sub.width = 5, sub.height = 5)
  pl.rank <- lapply(pl, function(x) return(x$p.rank))
  trash <- scFanc::wrap.plots.fanc(pl.rank, plot.out = paste0(out.dir, "/", root.name, "_allGroups_rank.png"), 
                          sub.width = 15, sub.height = 5, n.col = 1)
  return()
}
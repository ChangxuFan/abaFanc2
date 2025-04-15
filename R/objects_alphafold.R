#>>>>>>>> Generating Json from fasta.
alphafold.json.from.fa <- function(fa, ID.prefix = "", seeds = 10, json.dir) {
  seqs <- seqinr::read.fasta(fa, forceDNAtolower = F)
  ids <- alphafold.id.gen(length(seqs), prefix = ID.prefix)
  if (dir.exists(json.dir) || file.exists(json.dir)) {
    stop("dir.exists(json.dir) || file.exists(json.dir)")
  }
  dir.create(json.dir, showWarnings = F, recursive = T)
  seqs <- sapply(seqs, paste0, collapse = "")
  names(seqs) <- NULL
  
  if (length(seqs) != length(ids)) {
    stop("length(seqs) != length(ids)")
  }
  
  lapply(1:length(ids), function(i) {
    id <- ids[i]
    seq <- seqs[i]
    alphafold.jsongen(id = id, seq = seq, seeds = seeds, json.dir = json.dir)
  })
  return()
}

alphafold.jsongen <- function(id, seq, seeds = 10, json.dir = NULL) {
  jl <- list(
    name = id,
    modelSeeds = seeds,
    sequences = list(
      list(
        protein = list(
          id = id,
          sequence = seq
        )
      )
    ),
    dialect = "alphafold3",
    version = 2
  )
  
  json <- jsonlite::toJSON(jl, auto_unbox = T) %>% jsonlite::prettify() 
  if (!is.null(json.dir)) {
    write(json, paste0(json.dir, "/", id, '.json'))
  } else {
    cat(json)
  }
}

alphafold.id.gen <- function(n, prefix = "") {
  if (n > 26*26)
    stop("Currently only support n > 26")
  
  ids <- 1:n
  letter1 <- LETTERS[ceiling(ids/26)]
  letter2 <- ids %% 26
  letter2[letter2 == 0] <- 26
  letter2 <- LETTERS[letter2]
  
  out <- paste0(prefix, letter1, letter2) %>% toupper()
  if (any(duplicated(out))) {
    stop("any(duplicated(out))")
  }
  return(out)
}

#>>>>>>>>> parse json from MSA/template_search, and make the json for inference

alphafold.inference.json.gen <- function(jsons, out.json.dir, modelSeeds = 10, add.seed.to.name = F,
                                         rm.template = F, rm.MSA = F) {
  # these jsons must have one sequence per file.
  # if(dir.exists(out.json.dir) || file.exists(out.json.dir)) {
  #   if (overwrite) system(paste0("rm -rf ", out.json.dir))
  #   else stop("out.json.dir already exists")
  # }
  
  if (any(!file.exists(jsons))) {
    stop("any(!file.exists(jsons))")
  }
  
  sequences <- lapply(jsons, function(json) {
    jl <- jsonlite::read_json(json)
    if (length(jl$sequences) != 1) {
      warning(paste0("Json error: > 1 sequence in json. Taking the first one. \n", json))
    }
    
    seq <- jl$sequences[[1]]
    
    if (rm.MSA) {
      seq$protein$unpairedMsa <- ""
      seq$protein$pairedMsa <- ""
    }
    
    if (rm.template) {
      seq$protein$templates <- list()
    }
    return(seq)
  })
  
  names <- sapply(sequences, function(seq) seq[[1]]$id) %>% unname()
  name <- paste0(names, collapse = "_")
  
  if (length(modelSeeds) == 1) {
    modelSeeds <- list(modelSeeds)
  }
  
  jl <- list(
    dialect = "alphafold3",
    version = 2,
    name = name,
    sequences = sequences,
    modelSeeds = modelSeeds,
    bondedAtomPairs = NULL,
    userCCD = NULL
  )
  
  json <- jsonlite::toJSON(jl, auto_unbox = T, null = "null") %>% jsonlite::prettify() 
  
  dir.create(out.json.dir, showWarnings = F, recursive = T)
  
  if (length(modelSeeds) == 1 && add.seed.to.name)
    name <- paste0(name, "_seed", modelSeeds)
  write(json, paste0(out.json.dir, "/", name, '.json'))
  
  return()
}

alphafold.json.rm.slots <- function(json, out.dir = NULL, rm.MSA = F, rm.template = F,
                                    add.to.name = T) {
  if (rm.MSA + rm.template < 1) {
    stop("You aren't removing anything...")
  }
  jl <- jsonlite::read_json(json)
  
  jl$sequences <- lapply(jl$sequences, function(seq) {
    if (rm.MSA) {
      seq$protein$unpairedMsa <- ""
      seq$protein$pairedMsa <- ""
    }
    
    if (rm.template) {
      seq$protein$templates <- list()
    }
    return(seq)
  })
  
  if (is.null(out.dir)) out.dir <- dirname(json)
  root <- sub("(_data)*\\.json$", "", basename(json))
  dir.create(out.dir, showWarnings = F, recursive = T)
  suffix <- c("rmMSA", "rmTemplate")[c(rm.MSA, rm.template)] %>% 
    paste0(collapse = "_")
  out.json <- paste0(out.dir, "/", root, "_", suffix, ".json")
  if (add.to.name) jl$name <- paste0(jl$name, "_", suffix)
  
  json <- jsonlite::toJSON(jl, auto_unbox = T, null = "null") %>% jsonlite::prettify() 
  write(json, out.json)
  return(out.json)
}

# jl <- jsonlite::read_json("vl1_test_alphafold/res/vl1.5_alphafold3/rcsb_pdb_3d2u/rcsb_pdb_3d2u_rmB2M.json")

# alphafold.inference.json.gen("vl1_test_alphafold/res/vl1.5_alphafold3/b2m/b2m_data.json",
#                              out.json.dir = "vl1_test_alphafold/res/vl1.5_alphafold3/b2m/test_jsonOut/",
#                              overwrite = T)

alphafold.export.msa <- function(json) {
  jl <- jsonlite::fromJSON(json)
  out.dir <- paste0(sub(".json$", "", json), "_msa/")
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  lapply(1:nrow(jl$sequences$protein), function(i) {
    seq.df <- jl$sequences$protein[i, ]
    
    lapply(c("pairedMsa", "unpairedMsa"), function(slot) {
      msa <- seq.df[[slot]]
      if (is.null(msa)) return()
      
      outfile <- paste0(out.dir, "/", seq.df$id, "_", slot, ".fa")
      write(msa, outfile, sep = "\n")
    })
  })
  return()
}


alphafold.confidence.collect <- function(jsons, threads = 6, id.map.df = NULL,
                                         out.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(out.dir)
  
  exceptions <- which(!grepl("_summary_confidences.json", jsons))
  if (length(exceptions) > 0) {
    stop(paste0("Some of the jsons don't have the required suffix _summary_confidences.json: \n",
                paste0(jsons[exceptions[1:min(length(exceptions), 5)]])))
  }
  
  exceptions <- which(!file.exists(sub("_summary_confidences.json", "_data.json", jsons)))
  if (length(exceptions) > 0) {
    stop(paste0("Some of the jsons don't have the conresponding _data.json: \n",
                paste0(jsons[exceptions[1:min(length(exceptions), 5)]])))
  }
  
  jl <- utilsFanc::safelapply(jsons, function(json) {
    conf <- jsonlite::fromJSON(json)
    data <- jsonlite::fromJSON(sub("_summary_confidences.json", "_data.json", json))
    names <- data$sequences$protein$id
    if (!is.null(id.map.df)) {
      names <- alphafold.id.map.back(names, id.map.df = id.map.df)
    }
    
    name <- data$name
    
    conf$name <- name
    
    vec.slots <- c("chain_iptm", "chain_ptm")
    
    for (slot in vec.slots) {
      names(conf[[slot]]) <- names
    }
    
    mat.slots <- c("chain_pair_iptm", "chain_pair_pae_min")
    
    for (slot in mat.slots) {
      dimnames(conf[[slot]]) <- list(names, names)
      simple.mat <- conf[[slot]]
      simple.mat[lower.tri(simple.mat, diag = T)] <- NA
      
      melted <- reshape2::melt(
        simple.mat, value.name = stringr::str_extract(slot, "iptm|pae"), varnames = c("chain1", "chain2"))
      
      melted <- melted %>% na.omit()
      conf[[paste0(slot, ".melt")]] <- melted
    }

    return(conf)
  }, threads = threads)
  names(jl) <- sapply(jl, function(x) x$name) %>% unname()
  
  dir.create(out.dir, showWarnings = F, recursive = T)
  saveRDS(jl, paste0(out.dir, "/", root.name, "_jl.Rds"))
  
  vec.slots <- c("chain_iptm", "chain_ptm", "fraction_disordered", "has_clash",
                 "ptm", "ranking_score")
  
  utilsFanc::safelapply(vec.slots, function(slot) {
    df <- lapply(jl, function(conf) {
      df <- data.frame(V1 = conf[[slot]], V2 = conf$name)
      
      if (slot %in% c("chain_iptm", "chain_ptm")) {
        # vector slots
        names(df) <- c("chain", "file")
      } else {
        # scalar slots
        names(df) <- c(slot, "file")
      }
      return(df)
    }) %>% do.call(rbind, .)
    
    if (!is.null(id.map.df)) {
      df$int.name <- alphafold.id.map.back.int.name(df$file, id.map.df = id.map.df)
    }
    
    write.table(df, paste0(out.dir, "/", root.name, "_", slot, ".tsv"), quote = F,
                sep = "\t", row.names = F, col.names = T)
  }, threads = threads)
  
  mat.slots <- c("chain_pair_iptm.melt", "chain_pair_pae_min.melt")
  utilsFanc::safelapply(mat.slots, function(slot) {
    df <- lapply(jl, function(conf) conf[[slot]]) %>% do.call(rbind, .)
    df$file <- rownames(df) %>% sub("\\.\\d+", "", .)
    df <- utilsFanc::factor2character.fanc(df)
    
    if (!is.null(id.map.df)) {
      df$int.name <- alphafold.id.map.back.int.name(df$file, id.map.df = id.map.df)
    }
    
    write.table(df, paste0(out.dir, "/", root.name, "_", slot, ".tsv"), quote = F,
                sep = "\t", row.names = F, col.names = T)
  }, threads = min(threads, 2))
  
  return()
}

alphafold.confidence.mat.viz <- function(
    mat.file, out.suffix = NULL, id.map.df = NULL,
    hm.cutoff = 0.8,
    plot.col = "iptm", chain1.regex = NULL, chain2.regex = NULL) {
  
  df <- read.table(mat.file, header = T)
  required.cols <- c("chain1", "chain2", plot.col)
  utilsFanc::check.intersect(required.cols, "required.cols", colnames(df), "colnames(df)")
  
  if (!is.null(chain1.regex)) {
    df <- df %>% dplyr::filter(grepl(chain1.regex, chain1))
  }
  
  if (!is.null(chain2.regex)) {
    df <- df %>% dplyr::filter(grepl(chain2.regex, chain2))
  }
  
  if (nrow(df) < 1) {
    stop("After regex filtering there are 0 rows left")
  }
  
  utilsFanc::check.dups(paste0(df$chain1, "_", df$chain2), 'paste0(df$chain1, "_", df$chain2)')
  mat <- reshape2::acast(df, chain1 ~ chain2, value.var = plot.col)
  
  if (!is.null(id.map.df)) {
    colnames(mat) <- alphafold.id.map.back(colnames(mat), id.map.df = id.map.df)
    rownames(mat) <- alphafold.id.map.back(rownames(mat), id.map.df = id.map.df)
  }
  
  out.root <- sub(".tsv$", "", mat.file)
  if (!is.null(out.suffix)) out.root <- paste0(out.root, "_", out.suffix)
  
  saveRDS(mat, paste0(out.root, "_mat.Rds"))
  
  #------ plotting
  # heatmap:
  suppressMessages(library(ComplexHeatmap))
  p <- Heatmap(mat)
  print(p)
  scFanc::save.base.plot(p, file = paste0(out.root, "_hm.pdf"),
                         height = 20*nrow(mat) + 100, width = 20*ncol(mat) + 100)
  #rank plot
  
  if (!is.null(hm.cutoff)) {
    mat.cutoff <- mat
    if (plot.col == "iptm")
      mat.cutoff[mat.cutoff < hm.cutoff] <- 0
    if (plot.col == "pae")
      mat.cutoff[mat.cutoff > hm.cutoff] <- hm.cutoff
      
    p <- Heatmap(mat.cutoff)
    print(p)
    scFanc::save.base.plot(p, file = paste0(out.root, "_hm_", hm.cutoff, ".pdf"),
                           height = 20*nrow(mat) + 100, width = 20*ncol(mat) + 100)
  }
  
  lapply(c("chain1", "chain2"), function(chain.id) {
    if (chain.id == "chain1") mat <- t(mat)
    
    mat.df <- as.data.frame(mat)
    
    p <- scFanc::rank.plot(mat.df, vars = colnames(mat.df), na.last = FALSE, 
                           title = chain.id, use.lines = T, pt.size = 1)
    print(p)
    plot.out <- paste0(out.root, "_rank_", chain.id, ".pdf")
    ggsave(plot.out, p, device = cairo_pdf, width = 7, height = 5)
    
    # individual plots
    pl <- lapply(colnames(mat.df), function(var) {
      mat.df$names <- rownames(mat.df)
      p <- rank.plot(mat.df, vars = var, na.last = FALSE,
                             label.var = "names",
                             title = paste0(chain.id, " ", var),
                             use.lines = T, pt.size = 1)
    })
    p <- scFanc::wrap.plots.fanc(pl, plot.out = paste0(out.root, "_rank_", chain.id, "_each.pdf"))
    print(p)
    return()
  })
  
}


alphafold.id.map.back <- function(x, id.map.df) {
  # x: just a vector of alphafold names
  # id.map.df <- data.frame(alphafold = ..., real = ...)
  utilsFanc::check.intersect(c("alphafold", "real"), "required columns", 
                             colnames(id.map.df), "colnames(id.map.df)")
  
  lapply(id.map.df, function(x) {
    if (any(duplicated(x))) {
      stop("any(duplicated(x))")
    }
  })
  
  df <- data.frame(alphafold = x)
  df <- dplyr::left_join(df, id.map.df, by = "alphafold")
  df$real[is.na(df$real)] <- df$alphafold[is.na(df$real)]
  out <- df$real
  
  if (length(x) != length(out)) {
    stop("length(x) != length(out)")
  }
  return(out)
}

alphafold.id.map.back.int.name <- function(x, id.map.df) {
  # x <- KIRAA_HCMVAE_BM
  # x.df <- strsplit(x, "_") %>% as.data.frame() %>% t() %>% as.data.frame()
  # out <- lapply(x.df, alphafold.id.map.back, id.map.df = id.map.df) %>% as.list()
  # 
  # out <- Reduce(function(x, y) paste0(x, "_", y), out)
  out <- strsplit(x, "_") %>% 
    sapply(function(x) paste0(alphafold.id.map.back(x, id.map.df = id.map.df), collapse = "_")) %>% 
    unname()
  
  if (length(x) != length(out)) {
    stop("length(x) != length(out)")
  }
  return(out)
}

alphafold.multirun.comp <- function(
    mat.files, plot.col, hm.cutoff = 0.8,
    chain1.include = NULL, chain2.include = NULL,
    out.dir, root.name = NULL)
{
  if (is.null(root.name)) root.name <- basename(out.dir)
  out.root <- paste0(out.dir, "/", root.name, "_", plot.col)
  
  if (is.null(names(mat.files))) {
    stop("mat.files must be named")
  }
  
  df <- lapply(1:length(mat.files), function(i) {
    mat.file <- mat.files[i]
    
    df <- read.table(mat.file, header = T)
    required.cols <- c("chain1", "chain2", plot.col, "int.name")
    utilsFanc::check.intersect(required.cols, "required.cols", colnames(df), "colnames(df)")
    
    if (!is.null(chain1.include)) {
      df <- df %>% dplyr::filter(chain1 %in% chain1.include)
    }
    
    if (!is.null(chain2.include)) {
      df <- df %>% dplyr::filter(chain2 %in% chain2.include)
    }
    
    if (nrow(df) < 1) {
      stop("After inclusion filtering there are 0 rows left")
    }
    
    # utilsFanc::check.dups(paste0(df$chain1, "_", df$chain2), 'paste0(df$chain1, "_", df$chain2)')
    utilsFanc::check.dups(df$int.name, "df$int.name")
    df$run.name <- names(mat.files)[i]
    return(df)
  }) %>% do.call(rbind, .)
  
  mat <- reshape2::acast(df, int.name ~ run.name, value.var = plot.col)
  
  dir.create(dirname(out.root), showWarnings = F, recursive = T)
  saveRDS(mat, paste0(out.root, "_mat.Rds"))
  write.csv(mat, paste0(out.root, "_mat.csv"), quote = F)
  
  #------ plotting
  # heatmap:
  suppressMessages(library(ComplexHeatmap))
  p <- Heatmap(mat)
  print(p)
  scFanc::save.base.plot(p, file = paste0(out.root, "_hm.pdf"),
                         height = 20*nrow(mat) + 250, width = 20*ncol(mat) + 300)
  #rank plot
  
  if (!is.null(hm.cutoff)) {
    mat.cutoff <- mat
    if (plot.col == "iptm")
      mat.cutoff[mat.cutoff < hm.cutoff] <- 0
    if (plot.col == "pae")
      mat.cutoff[mat.cutoff > hm.cutoff] <- hm.cutoff
    
    p <- Heatmap(mat.cutoff)
    print(p)
    scFanc::save.base.plot(p, file = paste0(out.root, "_hm_", hm.cutoff, ".pdf"),
                           height = 20*nrow(mat) + 250, width = 20*ncol(mat) + 300)
  }
  return()
}

alphafold.plot.residue.mat <- function(residue.json, pae.cutoff = 2, contact.probs.cutoff = 0.8,
                                   out.dir, root.name = NULL,
                                   plot.height = 500, plot.width = 500) {
  
  # t <- jsonlite::read_json("vl5_multi_seeds/vl5.4.1_z_output/seed1/kirag_hcmvbu/kirag_hcmvbu_confidences.json", 
  #                          simplifyVector = T)
  # sapply(t, length)
  # # atom_chain_ids     atom_plddts   contact_probs             pae token_chain_ids   token_res_ids 
  # #           4273            4273             532             532             532             532 
  # 
  # tt <- jsonlite::fromJSON("vl5_multi_seeds/vl5.4.1_z_output/seed1/kirag_hcmvbu/kirag_hcmvbu_data.json")
  # tt$sequences$protein$sequence %>% nchar() %>% sum() # 532
  # 
  # # This means that each token is actually each amino acid
  # # token_res_ids seems to mean "token residue ids", 
  # # according to https://www.biorxiv.org/content/10.1101/2024.12.16.628070v1.full
  # 
  # mat <- t$contact_probs %>% do.call(rbind, .)
  # identical(mat, t(mat)) # TRUE
  
  if (is.null(root.name)) root.name <- basename(out.dir)
  out.root <- paste0(out.dir, "/", root.name, "_")
  
  jl <- jsonlite::fromJSON(residue.json)
  slots <- c("pae", "contact_probs")
  
  lapply(slots, function(slot) {
    mat <- jl[[slot]]
    
    lapply(c(F, T), function(bCutoff) {
      if (bCutoff) {
        cutoff <- ifelse(slot == "pae", pae.cutoff, contact.probs.cutoff)
        if (slot == "pae") mat[mat > cutoff] <- cutoff
        else mat[mat < cutoff] <- 0
      }
      
      chains <- unique(jl$token_chain_ids)
      anno.df <- data.frame(chain = jl$token_chain_ids)
      chain.colors <- utilsFanc::gg_color_hue(length(chains))
      names(chain.colors) <- chains
      anno.colors <- list(chain = chain.colors)
      
      left.annotation <- HeatmapAnnotation(df = anno.df, col = anno.colors,
                                           which = "row", show_legend = F, 
                                           annotation_name_rot = 0, show_annotation_name = F)
      top.annotation <- HeatmapAnnotation(df = anno.df, col = anno.colors,
                                          which = "column", show_legend = F,
                                          annotation_name_rot = 0, show_annotation_name = T)
      
      col_fun <- NULL
      if (slot == "pae") {
        col_fun = circlize::colorRamp2(c(min(mat), max(mat)), c("#03441B", "white"))
      }
      
      hm.params <- list(matrix = mat, name = slot, col = col_fun,
                        show_column_names = F, show_row_names = F,
                        show_row_dend = F, show_column_dend = F, 
                        cluster_columns = F, cluster_rows = F,
                        left_annotation = left.annotation, 
                        top_annotation = top.annotation,
                        width = unit(0.8, "snpc"), height = unit(0.8, "snpc"))
      hm <- do.call(what = ComplexHeatmap::Heatmap, args = hm.params)
      
      scFanc::save.base.plot(hm, file = paste0(out.root, slot, 
                                              ifelse(bCutoff, paste0("_cutoff", cutoff), ""),
                                              "_hm.png"),
                             height = plot.height, width = plot.width)
      
    })
    
  })
  
  
}






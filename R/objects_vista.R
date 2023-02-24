fa.2.vista <- function(fa, ref.id, query.id = NULL, out.file = NULL){
  # designed for pairend alignment results (2 sequences in fa, 1 ref 1 query)
  # if more sequences available, query.id must be specified.
  seqs <- seqinr::read.fasta(fa, forceDNAtolower = F)
  len <- sapply(seqs, length)
  if (length(unique(len)) > 1) {
    stop(paste0("sequences in ", fa, " have different lengths. Are they aligned?"))
  }
  if (! ref.id %in% names(seqs)) {
    stop("! ref.id %in% names(seqs)")
  }
  if (is.null(query.id) && length(seqs) == 2) {
    query.id <- names(seqs) %>% .[.!= ref.id]
  }
  if (! query.id %in% names(seqs)) {
    stop("! query.id %in% names(seqs)")
  }
  
  seqs <- seqs[c(ref.id, query.id)]
  dfs <- lapply(c(ref.id, query.id), function(name) {
    flag <- ifelse(name == ref.id, "r", "q")
    seq <- seqs[[name]]
    df <- data.frame(seq = as.character(seq), 
                     pos.align = 1:length(seq))
    df.nogap <- df %>% filter(seq != "-") %>% 
      dplyr::mutate(., pos.seq = 1:nrow(.))
    df <- left_join(df, df.nogap) %>% 
      dplyr::select(pos.align, pos.seq, seq)
    colnames(df) <- c("pos.align", paste0(flag, ".pos"), paste0(flag, ".seq"))
    return(df)
  })
  names(dfs) <- c("r", "q")
  df <- dplyr::left_join(dfs$r, dfs$q, by = "pos.align")
  df <- df %>% dplyr::mutate(
    dash = ifelse(r.seq == q.seq, "-", " "),
    r.seq = sub("\\-", "|", r.seq), q.seq = sub("\\-", "|", q.seq))
  df$r.pos[is.na(df$r.pos)] <- "   "
  df$q.pos[is.na(df$q.pos)] <- "   "
  
  out <- paste0(df$r.pos, "  ", df$r.seq, " ", df$dash, " ", df$q.seq, "  ", df$q.pos)
  if (is.null(out.file)) {
    out.file <- fa %>% tools::file_path_sans_ext() %>% paste0(".vistain")
  }
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  write(out, out.file, sep = "\n")
  return(out.file)
}

vista.run <- function(vistain.vec, ref.name, query.names = NULL, 
                      anno = NULL,
                      out.dir, root.name = NULL,
                      win.size = 100, resolution = 1,
                      start = NULL, end = NULL,
                      conservation.pid = 70, conservation.length = 100, write.align = T,
                      run = F,
                      vista.jar = "~/software/mVISTA/Vista.jar", 
                      pdf.jar = "~/software/mVISTA/retepPDF2.jar",
                      java = "/opt/apps/java/jdk1.8.0_92/bin/java") {
  # vistain.vec is a vector of alignment files in vista required format.
  # vistain.vec must be named. names(vistain.vec) will be used as sequence names in vista plot
  # added: if vistain.vec not named, you can supply names in query.names
  vistain.vec <- sapply(vistain.vec, normalizePath, mustWork = T)
  out.dir <- normalizePath(out.dir, mustWork = F)
  if (!is.null(query.names)) {
    if (length(vistain.vec) != length(query.names)) {
      stop("length(vistain.vec) != length(query.names)")
    }
    names(vistain.vec) <- query.names
  }
  if (is.null(names(vistain.vec))) {
    stop("vistain.vec must be named! you can also supply names via query.names")
  }
  
  if(is.null(root.name))
    root.name <- basename(out.dir)
  
  chunks <- paste0("TITLE", " ", root.name, "\n\n",
                   "OUTPUT", " ", out.dir, "/", root.name, ".pdf", "\n\n")
  for (q in names(vistain.vec)) {
    prefix <- paste0(out.dir, "/", root.name, "_", q, "_")
    chunks <- paste0(chunks, 
                     "ALIGN ", vistain.vec[q], "\n",
                     " SEQUENCES ", ref.name, " ", q, "\n",
                     " REGIONS ", conservation.pid, " ", conservation.length, "\n",
                     " REGION_FILE ", prefix, "cons_", conservation.pid, "_", conservation.length, "_regions.txt", "\n",
                     " SCORE_FILE ", prefix, "scores.txt", "\n")
    if (write.align) {
      chunks <- paste0(chunks, 
                       " ALIGNMENT_FILE ", prefix, "aln.txt", "\n")
    }
    chunks <- paste0(chunks, "END\n\n")
  }
  if (!is.null(anno)) {
    anno <- normalizePath(anno, mustWork = T)
    chunks <- paste0(chunks, "GENES ", anno, "\n\n")
  }
  chunks <- paste0(chunks,  "COORDINATE ", ref.name, "\n\n")
  if (!is.null(start)) {
    chunks <- paste0(chunks, "START ", start, "\n\n")
  }
  
  if (!is.null(end)) {
    chunks <- paste0(chunks, "END ", end, "\n\n")
  }
  
  chunks <- paste0(chunks, "RESOLUTION ", resolution, "\n\n",
         "WINDOW ", win.size, "\n\n",
         "NUM_PLOT_LINES ", 1, "\n\n",
         "LEGEND on \n\n",
         "AXIS_LABEL all \n\n")
  
  dir.create(out.dir, showWarnings = F, recursive = T)
  config.file <- paste0(out.dir, "/", root.name, ".config")
  write(chunks, config.file)
  cmd <- paste0("export CLASSPATH=", normalizePath(vista.jar, mustWork = T), ":",
                normalizePath(pdf.jar, mustWork = T),
                " && ", java, " Vista ", normalizePath(config.file, mustWork = T))
  if (run) {
    print(cmd); system(cmd)
  }
  cat(cmd)
  invisible(cmd)
}

vista.anno.write <- function(abao, bed, out.prefix = NULL) {
  # columns of bed: "chr", "start", "end", "name", "regionID"
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  if (is.null(out.prefix)) {
    out.prefix <- tools::file_path_sans_ext(bed) %>% paste0("_vistanno")
  }
  df <- utilsFanc::import.bed.fanc(bed = bed, no.shift = F)
  df <- df[, 1:5]
  colnames(df) <- c("chr", "start", "end", "name", "regionID")
  
  df <- df %>% dplyr::filter(regionID %in% abao@ori.df$regionID)
  if (nrow(df) < 1) {
    stop("nrow(df) < 1")
  }
  
  df %>% split(., f = df$regionID) %>% 
    lapply(function(df) {
      rID <- df$regionID[1]
      anno.1r <- df %>% split(., f = 1:nrow(df)) %>% 
        lapply(function(df) {
          ori <- abao@ori.df %>% dplyr::filter(regionID == rID)
          res <- data.frame(start = paste0("> ", df$start - ori$start + 1),
                            end = df$end - ori$start + 1,
                            name = df$name)
          return(res) 
        }) %>% Reduce(rbind, .)
      dir.create(dirname(out.prefix), showWarnings = F,recursive = T)
      out.file <- paste0(out.prefix, "_", rID, ".txt")
      write.table(anno.1r, out.file, sep = " ", quote = F, row.names = F, col.names = F)
      return(out.file)
    })
  return()
}

aba.read.slangan.core <- function(abao, seq1, seq2, SLangan.dir = "SLangan",
                                  write.vistain = F) {
  if (is.character(abao))
    abao <- readRDS(abao)
  if (!seq1 %in% abao@ori.df$regionID) {
    stop("!seq1 %in% abao@ori.df$regionID")
  }
  
  if (!seq2 %in% abao@ori.df$regionID) {
    stop("!seq2 %in% abao@ori.df$regionID")
  }
  
  aln.fa <- paste0(abao@work.dir, "/", SLangan.dir, "/", seq1, "_", seq2, ".fa")
  if (!file.exists(aln.fa)) {
    stop(paste0("alignment file does not exist: ", aln.fa))
  }
  seqs <- seqinr::read.fasta(aln.fa, 
                           forceDNAtolower = F)
  seqs <- seqs[c(seq1, seq2)]
  
  if (length(seqs[[1]])!=length(seqs[[2]])) {
    stop(paste0("sequence lengths differ: \n", 
                seq1, ": ", length(seqs[[seq1]]), "; ", seq2, ": ", length(seqs[[seq2]])))
  }
  ranges <- lapply(seqs, function(seq) {
    annot <- attr(seq, "Annot")
    if (!grepl("\\(\\+\\)", annot)) {
      stop("(+) was note detected in the header of sequence")
    }
    range <-  annot %>% sub(".+:", "", .) %>% sub(" .+$", "", .) %>% 
      strsplit("-") %>% unlist() %>% as.numeric()
    return(range)
  })
  rownames(abao@ori.df) <- abao@ori.df$regionID
  abao@ori.df$width <- abao@ori.df$end - abao@ori.df$start + 1
  widths <- abao@ori.df[c(seq1, seq2), "width"]
  names(widths) <- c(seq1, seq2)
  
  # now we pad the pairwise alignment:
  aln.df <- data.frame(seq1 = as.character(seqs[[seq1]]),
                       seq2 = as.character(seqs[[seq2]]))
  names(aln.df) <- c(seq1, seq2)
  
  pad.dfs <- lapply(1:2, function(i) {
    # i = 1: padded pseudo-alignment before real alignment starts. i = 2: after
    # we need to pad because shuffle langan doesn't always give you alignment that starts from 1.
    # this is how we pad it: if seq1 starts from 3, seq2 starts from 4:
    # seq1: NN---[real alignment]
    # seq2: --NNN[real alignment]
    # the same applies for the padding after alignment.
    if (i == 1) {
      pad.df <- data.frame(seq1 = c(rep("N", ranges[[seq1]][i] - 1), rep("-", ranges[[seq2]][i] - 1)))
    } else {
      pad.df <- data.frame(seq1 = c(rep("N", widths[seq1] - ranges[[seq1]][i]),
                                    rep("-", widths[seq2] - ranges[[seq2]][i])))
    }
    pad.df$seq2 <- rep(NA, length(pad.df$seq1))
    pad.df$seq2[pad.df$seq1 == "N"] <- "-"
    pad.df$seq2[pad.df$seq1 == "-"] <- "N"
    colnames(pad.df) <- c(seq1, seq2)
    return(pad.df)
  })
  
  aln.df <- rbind(pad.dfs[[1]], aln.df, pad.dfs[[2]])
  if (write.vistain) {
    vistain <- aln.fa %>% tools::file_path_sans_ext() %>% paste0(".vistain")
    seqinr::write.fasta(sequences = as.list(aln.df), names = colnames(aln.df), 
                        file.out = vistain)
    fa.2.vista(fa = vistain, out.file = vistain, ref.id = seq1)
  }
  
  invisible(aln.df)
}


vista.score.2.bdg <- function(score.file, chr, start, out.bdg = NULL, min.show = 0) {
  if (is.null(out.bdg)) {
    out.bdg <- tools::file_path_sans_ext(score.file) %>% paste0(".bedgraph")
  }
  
  scores <- read.table(score.file, sep = "\t")
  colnames(scores) <- c("start", "score")
  s <- start
  scores <- scores %>% 
    dplyr::mutate(chr = chr, start = start + s - 1, 
                  end = start, score = score - min.show) %>% 
    dplyr::filter(score > 0) %>% 
    dplyr::select(chr, start, end, score)
  utilsFanc::write.zip.fanc(df = scores, out.file = out.bdg, bed.shift = T)
  return(out.bdg)
}

vista.region.2.bed <- function(region.file, chr, start, out.bed = NULL,
                               must.exist = F) {
  if (is.null(out.bed)) {
    out.bed <- tools::file_path_sans_ext(region.file) %>% paste0(".bed")
  }
  if (!file.exists(region.file)) {
    if (must.exist) {
      stop(paste0("file not found: ", region.file))
    } else {
      df <- data.frame(chr = chr, start = start, end = start)
      df <- df[F, ]
    }
  } else {
    x <- readLines(region.file)
    x <- x[-c(1,2,3,4)]
    x <- x[x != ""]
    x <- x[!grepl("^Total", x)]
    x <- gsub(" +", "", x)
    s <- sub("\\(.+$", "", x) %>% as.numeric() 
    e <- sub(".+to", "", x) %>% sub("\\(.+$", "", .) %>% as.numeric() 
    df <- data.frame(chr = chr,
                     start = start + s - 1,
                     end = start + e -1)
  }
  utilsFanc::write.zip.fanc(df, out.file = out.bed, bed.shift = T)
  return(out.bed)
}
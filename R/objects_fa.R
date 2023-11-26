add.fa <- function(df, outdir) {
  system(paste0("mkdir -p ", outdir))
  df$fa <- paste0(outdir, "/", df$acc, ".fasta")
  for (i in 1:nrow(df)) {

    fa <- entrez_fetch("nuccore", id = df[i, "acc"], rettype = "fasta", api_key ="4d2374d16cee57dded1296ac1a48ba9c3b09")
    write(fa, df[i, "fa"])
  }

  trash <- lapply(df$fa, function(x) {
    if(!file.exists(x)) {
      stop(paste0(x, " was not successfully created"))
    }
  })
  return(df)
}

get.fasta.bed <- function(bed, root.name=NULL, genome=NULL, source.fasta=NULL, add.coordinate=T, add.additional.columns=T, return.fasta = F) {

  if (is.null(root.name)) {
    if (is.character(bed))
      root.name <- sub(".bed", "", bed)
    else
      stop("root.name not specified")
  }

  if (is.character(bed))
    bed <- read.table(bed, as.is = T, sep = "\t")

  bed.out <- bed %>% mutate(id=1:nrow(.)) %>% split(., f=.$id) %>% lapply(function(x) {
    x$id <- NULL
    name.col <- ""
    if (add.coordinate == T)
      name.col <- paste0(x[1,1], ":", x[1,2], "-", x[1,3])
    if (ncol(x) >= 4 && add.additional.columns==T) {
      if (name.col != "")
        name.col <- paste0(name.col, "|")
      name.col <- name.col %>% paste0(paste0(x[1,4:ncol(x)], collapse = "|"))
    }
    x[1,4] <- name.col
    return(x[,1:4])
  }) %>% Reduce(rbind,.)

  write.table(bed.out, paste0(root.name, ".pre.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
  out.fa <- paste0(root.name, ".fa")
  if (!is.null(genome)) {
    fi <- paste0("/bar/cfan/genomes/", genome, "/", genome, ".fa")
  }
  if (!is.null(source.fasta)) {
    fi <- source.fasta
  }
  cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/bedtools getfasta -fi ",fi,
    " -fo ", out.fa, " -bed ", root.name, ".pre.bed -name")
  print(cmd)
  system(cmd)

  fasta.wrap.fanc(in.fa = out.fa)
  if (return.fasta == T) {
    return(seqinr::read.fasta(out.fa))
  }
  return(out.fa)

}

get.fasta.bed.2 <- function(df, genome=NULL, fa=NULL,  df.is.bed, threads = 1, out.fa = NULL) {
  # df should contain: chr, start, end, regionID, fa and/or genome. strand is optional
  # if all regions share the same genome or fa, they can also be specified through arguments
  if (is.character(df)) {
    if (df.is.bed == T) {
      df.bed <- read.table(df, sep = "\t", as.is = T)
      df <- df.bed[, 1:4]
      colnames(df) <- c("chr", "start", "end", "regionID")
      if (!is.null(df.bed$V6))
        df$strand <- df.bed$V6
    } else {
      df <- read.table(df, sep = "\t", as.is = T, quote = "", header = T)
    }
  }

  if (!is.null(genome))
    df$genome <- genome
  if (!is.null(fa))
    df$fa <- fa

  if (sum(duplicated(df$regionID)) > 0)
    stop("regionID must be unique")

  seq.list <- df %>% split(., f = factor(.$regionID, levels = .$regionID)) %>%
    mclapply(function(x) {
      if (is.null(x$fa) || is.na(x$fa)) {
        x$fa <- paste0("/bar/cfan/genomes/", x$genome, "/", x$genome, ".fa")
      }
      seq <- get.fasta.core(genome.fa = x$fa, chr = x$chr, start = x$start,
                            end = x$end, strand = x$strand)
      return(seq)
    }, mc.cores = threads, mc.cleanup = T)
  if (!is.null(out.fa)) {
    seqinr::write.fasta(sequences = seq.list, names = names(seq.list), file.out = out.fa)
  }
  return(seq.list)
}

get.fasta.core <- function(genome.fa, chr, start, end, strand=NULL,
                           bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  df <- data.frame(chr = chr, start = start, end = end)
  if (!is.null(strand))
    df$strand <- strand
  if (nrow(df) != 1)
    stop("get.fasta.core: df must be exactly 1 row")
  gr <- GenomicRanges::makeGRangesFromDataFrame(df)
  bed <- tempfile()
  rtracklayer::export.bed(gr, bed)
  fa <- tempfile()
  cmd <- paste0(bedtools, " getfasta -fi ", genome.fa, " -bed ", bed,
                " -fo ", fa)
  if (!is.null(strand))
    cmd <- paste0(cmd, " -s")
  system(cmd)
  fa <- seqinr::read.fasta(fa, forceDNAtolower = F)
  seq <- seqinr::getSequence(fa)[[1]]
  return(seq)
}


get.fasta.gr <- function(gr, id.col = NULL, genome, fa = NULL, outfile = NULL,
                          drop.strand = F,
                         wrap.fa = F, return.fa = T, as.string = T, print.cmd = T,
                         bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  # written: 2022-01-10. latest function of this get fasta family.
  # note I validated this function through eyeballing chr6:130124429-130139101 to make sure the 2 ends are exact
  if (is.null(fa))
    fa <- paste0("~/genomes/", genome, "/", genome, ".fa")
  if (!file.exists(fa)) {
    stop("!file.exists(fa)")
  }
  if (is.null(id.col)) {
    id.col <- "id"
    gr$id <- paste0("region_", 1:length(gr), "|",utilsFanc::gr.get.loci(gr), strand(gr))
  }
  id.vec <- mcols(gr)[, id.col] %>% as.character()
  mcols(gr) <- NULL
  gr$id <- id.vec

  bed <- tempfile()
  if (is.null(outfile))
    outfile <- tempfile()
  dir.create(dirname(outfile), showWarnings = F, recursive = T)
  trash <- utilsFanc::write.zip.fanc(df = gr, out.file = bed, bed.shift = T, zip = F)
  s <- "-s"
  if (drop.strand)
    s <- ""
  cmd <- paste0(bedtools, " getfasta ", s, " -fi ", fa, " -fo ", outfile, " -bed ", bed, " -name")
  if (print.cmd == T)
    print(cmd)
  system(cmd)
  if (wrap.fa == T)
    abaFanc2::fasta.wrap.fanc(in.fa = outfile)
  if (return.fa == T) {
    res.fa <- seqinr::read.fasta(outfile, as.string = as.string, forceDNAtolower = F)
    names(res.fa) <- names(res.fa) %>% sub("\\([\\+\\-]\\)$", "", .)
    if (! identical(names(res.fa), id.vec )) {
      stop("! identical(names(res.fa), id.vec )")
    }
    return(res.fa)
  }
  return()
}

label.gaps <- function(fa, out.bed = NULL) {
  if (is.null(out.bed)) {
    out.bed <- paste0(tools::file_path_sans_ext(fa), "_gaps.bed")
  }
  seqs <- seqinr::read.fasta(fa)

  gr <- lapply(seqs, function(seq) {
    df <- data.frame(chr = attr(seq,"name"),
                     start = 1:length(seq), end = 1:length(seq),
                     seq = as.vector(seq) %>% toupper())
    df <- df %>% dplyr::filter(seq == "N") %>% dplyr::mutate(seq = NULL)
    gr <- makeGRangesFromDataFrame(df) %>% GenomicRanges::reduce()
    return(gr)
  }) %>% Reduce(c, .)
  gr$type <- "gap"
  utilsFanc::write.zip.fanc(gr, out.bed, bed.shift = T)
  return(out.bed)
}

fa.mask.regions <- function(fa.vec, gr.list, master.dir, bowtie2_index = T) {
  # for each gr in gr.list, all fa's in fa.vec will be masked
  # gr.list must be named
  if (is.null(names(gr.list))) {
    stop("is.null(names(gr.list))")
  }
  lapply(names(gr.list), function(region) {
    gr <- gr.list[[region]]
    lapply(fa.vec, function(fa) {
      seqs <- seqinr::read.fasta(fa, forceDNAtolower = F)
      seqs <- lapply(seqs, function(seq) {
        for (i in 1:length(gr)) {
          line.gr <- gr[i]
          if (as.character(seqnames(line.gr)) == attr(seq, "name")) {
            seq[start(line.gr):end(line.gr)] <- "N"
          }
        }
        return(seq)
      })
      root.name <- paste0(basename(fa) %>% tools::file_path_sans_ext(), "_", region)
      out.fa <- paste0(master.dir, "/", root.name, "/", root.name, ".fa")
      system(paste0("mkdir -p ", dirname(out.fa)))
      seqinr::write.fasta(sequences = seqs, names = names(seqs),
                          file.out = out.fa)
      label.gaps(fa = out.fa)
      system(paste0(SAMTOOLS, " faidx ", out.fa))
      if (bowtie2_index == T) {
        liteRnaSeqFanc::bowtie2.index(out.fa, threads = 1)
      }
      return()
    })
    return()
  })
  return()
}

fa.remove.repeat <- function(in.fa) {
  fa <- seqinr::read.fasta(file = in.fa, forceDNAtolower = F)
  map.dir <- paste0(in.fa, "_map/")
  dir.create(map.dir, showWarnings = F, recursive = T)
  noR.seqs <- lapply(fa, function(seq) {
    name <- attr(seq, "name")
    seq <- as.character(seq)
    map <- data.frame(ori.seq = seq, ori.pos = 1:length(seq))
    map <- map %>% dplyr::filter(ori.seq %in% LETTERS)
    if (nrow(map) < 1) {
      stop("nrow(map) < 1")
    }
    map$noRepeat.pos <- 1:nrow(map)
    write.table(map, file = paste0(map.dir, "/", name, "_map.tsv"), row.names = F,
                col.names = T, quote = F, sep = "\t")
    return(map$ori.seq)
  })
  names(noR.seqs) <- names(fa)
  out.fa <- utilsFanc::insert.name.before.ext(name = in.fa,insert = "noR", delim = "_")
  seqinr::write.fasta(sequences = noR.seqs, names = names(noR.seqs), file.out = out.fa)
  return(out.fa)
}

fa.mask.repeat <- function(in.fa, out.fa = NULL) {
  if (is.null(out.fa)) {
    out.fa <- utilsFanc::insert.name.before.ext(name = in.fa, insert = "maskR", delim = "_")
  }
  fa <- seqinr::read.fasta(file = in.fa, forceDNAtolower = F)
  fa <- lapply(fa, function(fa) {
    bk <- fa
    fa[fa %in% base::letters] <- "n"
    return(fa)
  })
  dir.create(dirname(out.fa), showWarnings = F, recursive = T)
  seqinr::write.fasta(sequences = fa, names = names(fa), file.out = out.fa)
}
fa.aln.remove.gap <- function(in.fa, out.fa = NULL) {
  fa.mat <- Biostrings::readDNAMultipleAlignment(in.fa) %>% as.matrix()
  b.gap.mat <- (!fa.mat %in% c("A", "T", "C", "G", "a", "t", "c", "g")) %>% matrix(ncol = ncol(fa.mat))
  b.keep <- colSums(b.gap.mat) == 0
  if (is.null(out.fa)) {
    out.fa <- utilsFanc::insert.name.before.ext(name = in.fa, insert = "noGap", delim = "_")
  }
  dir.create(path = dirname(out.fa), showWarnings = F, recursive = T)

  map <- data.frame(ori.pos = 1:ncol(fa.mat), b.keep = b.keep) %>% filter(b.keep == T) %>%
    dplyr::mutate(., new.pos = 1:nrow(.), b.keep = NULL)
  write.table(map, paste0(out.fa, "_map.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  new.fa.mat <- fa.mat[, b.keep]

  ape::write.dna(x = new.fa.mat, file = out.fa, format = "fasta", colsep = "")
  fasta2phylip(out.fa)
  return(out.fa)
}

fa.rm.dot <- function(fa) {
  # remove the dot in sequences names
  fasta <- seqinr::read.fasta(fa)
  names(fasta) <- names(fasta) %>% gsub("\\.", "_", .)
  seqinr::write.fasta(sequences = fasta, names = names(fasta), file.out = fa)
  if (file.exists(sub("\\.fa(sta)*", ".phy", fa))) {
    fasta2phylip(fa)
  }
  return()
}

fa.2.nex <- function(fa.vec, out.dir = NULL, out.file = NULL, 
                     use.partition = T, add.mb.GTR = F, ngen = 100000) {
  # partition model for MrBayes
  # fa.vec must be named
  if (is.null(names(fa.vec)) && use.partition) {
    stop("is.null(names(fa.vec))")
  }
  mat.list <- lapply(fa.vec, function(fa) {
    mat <- Biostrings::readDNAMultipleAlignment(fa) %>% Biostrings::as.matrix()
    return(mat)
  })
  taxons <- lapply(mat.list, rownames) %>% Reduce(union,. )
  mat.list <- lapply(mat.list, function(mat) {
    add <- taxons %>% .[!.%in% rownames(mat)]
    if (length(add) > 0) {
      add.mat <- matrix("?", ncol = ncol(mat), nrow = length(add))
      rownames(add.mat) <- add
      mat <- rbind(mat, add.mat)
    }
    mat <- mat[taxons, ]
    return(mat)
  })
  len <- sapply(mat.list, ncol)
  mat <- Reduce(cbind, mat.list)
  if (is.null(out.file) && length(fa.vec) == 1) {
    if (is.null(out.dir)) {
      out.dir <- dirname(fa.vec) %>% paste0("/mb/")
    }
    root <- basename(fa.vec) %>% tools::file_path_sans_ext()
    out.file <- paste0(out.dir, "/", root, ".nex")
  }
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  ape::write.nexus.data(ape::as.DNAbin(mat), file = out.file)
  if (use.partition) {
    write("begin mrbayes;", file = out.file, append = T)
    total <- 0
    for (i in 1:length(len)) {
      start <- total + 1
      end <- len[i] + start - 1
      cmd <- paste0("    charset ", names(len)[i], " = ", start, "-", end, ";")
      write(cmd, file = out.file, append = T)
      total <- total + len[i]
    }
    write(paste0("    partition favored = ", length(len), ": ",
                 paste0(names(len), collapse = ", "), ";"), file = out.file, append = T)
    write("    set partition=favored;", file = out.file, append = T)
    if (add.mb.GTR) {
      cmd <- c(paste0("lset app=(", paste0(1:length(len), collapse = ","), ") rates=invgamma nst=6;"),
               "unlink revmat=(all) pinvar=(all) shape=(all) statefreq=(all);",
               "prset applyto=(all) ratepr=variable;",
               paste0("mcmc ngen=", ngen,";"),
               "sump;",
               "sumt;") %>% paste0("    ", .)
      write(cmd, file = out.file, append = T, sep = "\n")
    }

    write("end; ", file = out.file, append = T)

  } else if (add.mb.GTR) {
    write("begin mrbayes;", file = out.file, append = T)
    cmd <- c("lset nst=6 rates=invgamma;",
             paste0("mcmc ngen=", ngen,";"),
             "sump;",
             "sumt;") %>% paste0("    ", .)
    write(cmd, file = out.file, append = T, sep = "\n")
    write("end; ", file = out.file, append = T)
  }
  
  return(out.file)
}

nex.2.fa <- function(nex, out.fa = NULL) {
  if (is.null(out.fa)) {
    out.fa <- paste0(sub(".nex$", "", nex), ".fa")
  }
  d <- ape::read.nexus.data(file = nex)
  dir.create(dirname(out.fa), showWarnings = F, recursive = T)
  seqinr::write.fasta(sequences = d, names = names(d),
                      file.out = out.fa)
  return()
}

fa.simu.sv <- function(in.fa, sv.bed) {
  # format for sv.bed
  # chr, start (0 based), end, sv.type, additional info
  fa <- seqinr::read.fasta(in.fa) %>% lapply(as.vector)
  ids <- lapply(fa, function(x) return(1:length(x)))
  sizes <- lapply(ids, function(x) return(x[length(x)]))
  if (is.character(sv.bed)) {
    sv.df <- utilsFanc::import.bed.fanc(bed = sv.bed)
  } else {
    sv.df <- sv.bed
  }
  sv.df <- sv.df[, 1:5]
  colnames(sv.df) <- c("chr", "start", "end", "sv.type", "sv.info")
  utilsFanc::check.intersect(sv.df$chr, "sv.df$chr", names(fa), "names(fa)")

  for (i in 1:nrow(sv.df)) {
    sv <- sv.df[i, ]
    id <- ids[[sv$chr]]
    if (sv$start > sv$end + 1) {
      stop(paste0("error at ", i, ": sv$start > sv$end + 1"))
    }
    if (sv$end > max(id)) {
      stop(paste0("error at ", i, ": sv$end > max(id)"))
    }
    if (sv.type == "deletion") {
      to.del <- sv$start:sv$end
      id <- ! id %in% to.del
    }
    if (sv.type == "insertion") {
      if (sv$start - sv$end < 1) {
        # this is how you code an insertion in bed: chr start start, after bed shift this will be
        # chr start+1 start
        # so a pure insertion will have sv$start - sv$end == 1, but you could also remove
        # some sequences at the insertion junction, which leads to sv$start < sv$end
        to.del <- sv$start:sv$end
        id <- ! id %in% to.del
      }
      if (grepl(":\\d+\\-\\d+", sv$sv.info)) {
        # this means that inserted sequence is encoded as a locus. grep this sequence from fa
        # multiple loci can be supplied via ";" deliminator.
        loci <- sv$sv.info %>% strsplit(";") %>% unlist() %>%
          utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = T)
        if (length(unique(loci$chr)) != 1)
          stop("length(unique(loci$chr)) != 1")
        if (loci$chr[1] != sv$chr) {
          stop("loci$chr[1] != sv$chr")
        }
        for (j in 1:nrow(loci)) {
          to.ins <- loci$start[j]:loci$end[j]
          # it's not entirely easy to find where to insert:
          # sv$start and sv$end might have been duplicated or inverted or whatever...
          # but we look for the closest 2 so that we remove
        }
      } else if (grepl("^[ATCGatcgNn]+$", sv$sv.info)) {
        stop("not developed yet")
      } else {
        stop(paste0("error at ", i, ": insertion's sv.info is not of correct format"))
      }
    }

    ids[[sv$chr]] <- id
  }

}

fa.sep <- function(in.fa, out.prefix = NULL) {
  if (is.null(out.prefix)) {
    out.prefix <- in.fa %>% tools::file_path_sans_ext()
    out.prefix <- paste0(out.prefix, "/", basename(out.prefix))
  }
  seqs <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
  dir.create(dirname(out.prefix), showWarnings = F, recursive = T)
  lapply(names(seqs), function(name) {
    seq <- seqs[[name]]
    seqinr::write.fasta(list(seq), names = list(name),
                        file.out = paste0(out.prefix, "_", name, ".fa"))
  })
  return()
}

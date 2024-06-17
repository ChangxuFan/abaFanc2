# #
SHORT.ALIGN.PARAMS <- " -word_size 7 -gapopen 5 -gapextend 2 -max_target_seqs 5000 -evalue 10"

fasta.wrap.fanc <- function(in.fa, out.fa=NULL) {
  if (is.null(out.fa))
    out.fa <- in.fa
  fa <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
  headers <- lapply(fa, attr, "Annot") %>% unlist() %>% sub("^>", "", .)
  seqinr::write.fasta(fa, headers, out.fa)
  return(out.fa)
}

read.fasta.fanc <- function(in.fa, return.list = T) {
  fa <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
  if (return.list == T)
    fa <- seqinr2list(fa)
  return(fa)
}

seqinr2list <- function(seqinr.o) {
  names <- names(seqinr.o)
  fa.list <- seqinr::getSequence(seqinr.o)
  names(fa.list) <- names
  return(fa.list)
}
# fasta.mask.fanc <- function(in.fa, out.fa=NULL) {
#   if (is.null(out.fa))
#     out.fa <- in.fa
#   fa <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
#   headers <- lapply(fa, attr, "Annot") %>% unlist() %>% sub("^>", "", .)
#   seqs <- lapply(fa, function(x) paste0(x[1:length(x)], collapse = "") ) %>%
#     unlist() %>% `names<-`(NULL) %>%
#     gsub("[a-z]", "", .)
#   seqinr::write.fasta(strsplit(seqs, ""), headers, out.fa)
#   return(out.fa)
# }


# bash2ftp <- function(filename) {
#   ftp <- sub("^~", "https://wangftp.wustl.edu/~cfan", filename) %>%
#     sub("/bar/cfan", "https://wangftp.wustl.edu/~cfan", .)
#   return(ftp)
# }

# bed.align <- function(bed, root.name=NULL, genome, align=F, mafft.options = "--auto --reorder",
#  add.coordinate=T, mask =F) {

#   if (is.null(root.name)) {
#     if (is.character(bed))
#       root.name <- sub(".bed", "", bed)
#     else
#       stop("root.name not specified")
#   }

#   if (is.character(bed))
#     bed <- read.table(bed, as.is = T, sep = "\t")

#   bed.out <- bed %>% mutate(id=1:nrow(.)) %>% split(., f=.$id) %>% lapply(function(x) {
#     x$id <- NULL
#     name.col <- ""
#     if (add.coordinate == T)
#       name.col <- paste0(x[1,1], ":", x[1,2], "-", x[1,3])
#     if (ncol(x) >= 4) {
#       if (name.col != "")
#         name.col <- paste0(name.col, "|")
#       name.col <- name.col %>% paste0(paste0(x[1,4:ncol(x)], collapse = "|"))
#     }
#     x[1,4] <- name.col
#     return(x[,1:4])
#   }) %>% Reduce(rbind,.)
  
#   write.table(bed.out, paste0(root.name, ".pre.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
#   cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/bedtools getfasta -fi /bar/cfan/genomes/",
#                 genome, "/", genome, ".fa -fo ", root.name, ".pre.fa -bed ", root.name, ".pre.bed -name")
#   print(cmd)
#   system(cmd)
  
#   fasta.wrap.fanc(in.fa = paste0(root.name, ".pre.fa"))

#   to.mafft <- paste0(root.name, ".pre.fa")

#   if (mask == T) {
#     masked <- fasta.mask.fanc(in.fa = to.mafft, out.fa = paste0(to.mafft, ".masked"))
#     to.mafft <- paste0(to.mafft, ".masked")
#   } 

#   if (align == T) {
#     after.mafft <- paste0(root.name, ".aligned.fa")
#     cmd <- paste0("/opt/apps/mafft/7.427/mafft ", mafft.options, " ", to.mafft, " > ", after.mafft)
#     print(cmd)
#     system(cmd)
#   } else {
#     after.mafft <- to.mafft
#   }

#   return(after.mafft)

# }


mafft.fanc <- function(in.fa, tempt.fa = NULL, mafft = "/opt/apps/mafft/7.427/mafft", aligned.fa = NULL,
                       mafft.params = "--auto --reorder --distout --treeout") {
  if (!is.character(in.fa)) {
    if (is.null(tempt.fa))
      in.fa.file <- tempfile() 
    else
      in.fa.file <- tempt.fa
    seqinr::write.fasta(in.fa, names(in.fa), in.fa.file)
    in.fa <- in.fa.file
  }
  if (is.null(aligned.fa))
    aligned.fa <- tempfile()
  cmd <- paste0(mafft, " ", mafft.params, " ", in.fa,  " > ", aligned.fa)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
  out.fa <- seqinr::read.fasta(aligned.fa, forceDNAtolower = F)
  names <- names(out.fa)
  out.fa <- seqinr::getSequence(out.fa)
  names(out.fa) <- names
  return(out.fa)
}

mafft.fanc.pairwise <- function(in.fa, ref, 
                                tempt.dir = NULL, mafft = "/opt/apps/mafft/7.427/mafft", 
                                aligned.fa = NULL, mafft.params = "--auto", threads = 1,
                                double.check = F) {
  seqinr.in <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
  in.seqs <- seqinr::getSequence(seqinr.in)
  names(in.seqs) <- seqinr::getName(seqinr.in)
  rm(seqinr.in)
  
  if (any(duplicated(names(in.seqs)))) {
    stop("duplicated seq names found in in.fa")
  }
  
  if (! ref %in% names(in.seqs)) {
    stop("ref is not found in in.fa")
  }
  
  ref.seq <- in.seqs[[ref]]
  ref.length <- length(ref.seq)
  if (is.null(tempt.dir)) tempt.dir <- tempdir()
  tempt.dir <- paste0(tempt.dir, "/ref_", ref, "/")
  
  dir.create(paste0(tempt.dir), showWarnings = F, recursive = T)
  
  queries <- names(in.seqs) %>% .[. != ref]
  aligned.raw <- utilsFanc::safelapply(queries, function(query.id) {
    in.fa <- paste0(tempt.dir, "/", query.id, ".fa")
    aligned.fa <- paste0(tempt.dir, "/", query.id, "_aln.fa")
    seqinr::write.fasta(sequences = in.seqs[c(ref, query.id)], names = c(ref, query.id),
                        file.out = in.fa)
    cmd <- paste0(mafft, " ", mafft.params, " ", in.fa,  " > ", aligned.fa)
    utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
    out.fa <- seqinr::read.fasta(aligned.fa, forceDNAtolower = F)
    names <- names(out.fa)
    out.fa <- seqinr::getSequence(out.fa)
    names(out.fa) <- names
    out.fa <- out.fa[c(ref, query.id)]
    return(out.fa)
  }, threads = threads)
  
  aligned.gapLabel <- lapply(aligned.raw, function(df) {
    df <- as.data.frame(df)
    df$pos.aln <- 1:nrow(df)
    df <- rbind(df[df[, ref] != "-",], df[df[, ref] == "-",])
    df$pos.ref <- 1:nrow(df)
    df$pos.ref[df[, ref] == "-"] <- 0
    if (max(df$pos.ref) != ref.length) {
      stop("max(df$pos.ref) != ref.length")
    }
    df <- df %>% dplyr::arrange(pos.aln)
    x <- rle(df$pos.ref)
    gaps <- which(x$values == 0)
    if (gaps[1] == 1) {
      gaps <- gaps[-1]
    }
    x$values[gaps] <- x$values[gaps - 1]
    df$pos.ref <- inverse.rle(x)
    df$pos.ref <- make.unique(as.character(df$pos.ref))
    df$pos.aln <- NULL
    return(df)
  })
  
  aligned <- Reduce(function(x, y) dplyr::full_join(x, y, by = c(ref, "pos.ref")), aligned.gapLabel)
  
  aligned$pos.int <- sub("\\..+$", "", aligned$pos.ref) %>% as.numeric()
  aligned$pos.dec <- sub("^.+\\.", "", aligned$pos.ref) %>% as.numeric()
  aligned$pos.dec[!grepl("\\.", aligned$pos.ref)] <- 0
  
  aligned <- aligned %>% arrange(pos.int, pos.dec)
  
  aligned[, c("pos.int", "pos.dec", "pos.ref")] <- NULL
  aligned <- as.matrix(aligned)
  
  aligned[is.na(aligned)] <- "-"
  aligned <- aligned[, c(ref, queries)]
  if (double.check) {
    # this part is to make sure that everything is correct: 
    lapply(names(in.seqs), function(seq.name) {
      ori <- in.seqs[[seq.name]] %>% tolower()
      final <- aligned[, seq.name] %>% .[.!="-"] %>% tolower()
      if (!identical(final, ori)) {
        stop(paste0("sequence doesn't match for: ", seq.name))
      }
    })
    
    # make sure that the pairwise alignments are also preserved:
    lapply(aligned.raw, function(aln.raw) {
      aln.raw <- aln.raw %>% as.data.frame()
      aln.final <- aligned[, c(names(aln.raw))]
      aln.final <- aln.final[!(aln.final[, 1] == "-" & aln.final[, 2] == "-"), ] %>% 
        as.data.frame()
      
      aln.final <- aln.final %>% lapply(tolower) %>% as.data.frame()
      aln.raw <- aln.raw %>% lapply(tolower) %>% as.data.frame()
      # print(head(aln.raw))
      # print(head(aln.final))
      if (!identical(aln.final, aln.raw)) {
        stop(paste0("pairwise alignment couldn't be reproduced for: ", 
                    paste0(names(aln.raw), sep = ", ")))
      }
    })
  }
  aligned <- as.list(as.data.frame(aligned))
  
  if (!is.null(aligned.fa)) {
    dir.create(dirname(aligned.fa), showWarnings = F, recursive = T)
    seqinr::write.fasta(aligned, names = names(aligned),
                        file.out = aligned.fa)
  }
  return(aligned)
}



gap.summary <- function(x) {
  # x needs to be a vector
  # try: gap.summary(unlist(strsplit("---ATCG--A-", "")))
  x <- x == "-"
  r <- rle(x)
  non.gap.length <- r$lengths
  non.gap.length[r$values] <- 0
  non.gap.cumsum <- cumsum(non.gap.length)
  summ <- data.frame(gap.pos = non.gap.cumsum[r$values], gap.length = r$lengths[r$values])
  return(summ)
}


prank.fanc <- function(in.fa, out.root.name = NULL, add.params = "", use.F = T, 
                       prank = PRANK, log.file = NULL) {
  if (is.null(out.root.name)) {
    out.root.name <- paste0(dirname(in.fa), "/prank/",
                            basename(in.fa) %>% sub("(in)*.fa", "prank",.))
    if (use.F) {
      out.root.name <- paste0(out.root.name, "_F")
    } else {
      out.root.name <- paste0(out.root.name, "_Fless")
    }
  }
  dir.create(dirname(out.root.name), showWarnings = F, recursive = T)
  cmd <- paste0(prank, " -d=", in.fa, " -o=", out.root.name, " -showanc -showiter")
  if (use.F == T)
    cmd <- paste0(cmd, " -F")
  if (is.null(log.file))
    log.file <- paste0(out.root.name, ".log")
  utilsFanc::cmd.exec.fanc(cmd, run = T, stdout.file = log.file, intern = F)
  if(length(Sys.glob(out.root.name %>% paste0("*"))) < 1)
    stop(paste0(out.root.name, " was not successfully generated"))
  return(paste0(out.root.name, ".best.fas"))
}

# mpw <- function(from.fa=NULL, from.bed=NULL, ref.pattern, chr.name, shift=0,
#                 genome=NULL, tempt.dir, SNV.dir, aligner, dry=T, out.json=F) {
#   # note: shift can be set at equal to the left of the reference region in the bed file.
#   if (is.null(from.fa) && is.null(from.bed))
#     stop("at least one of fa or bed should be supplied")
#   if (!is.null(from.fa) && !is.null(from.bed))
#     stop("only one of fa and bed can be supplied, you supplied both somehow")
#   system(paste0("rm -rf ", tempt.dir, " ", SNV.dir))
#   system(paste0("mkdir -p ", tempt.dir, " ", SNV.dir))
#   if (!is.null(from.bed)) {
#     if (is.null(genome))
#       stop("if bed is supplied, genome must also be specified")
#     from.fa <- bed.align(bed = from.bed,  root.name = paste0(tempt.dir, "/mpw"), genome = genome, align = F)
#   }

#   # now extract reference from the fasta file:
#   fa <- seqinr::read.fasta(from.fa)
#   fa.ref <- grep(ref.pattern,sapply(fa, attr, "Annot"))
#   if (length(fa.ref) > 1)
#     stop("ref.pattern matched to more than one sequence, they are: " %>% paste0(paste0(fa.ref, collapse = ", ")))
#   if (length(fa.ref) == 0)
#     stop("ref.pattern didn't match any thing")

#   fa.ref <- fa[[fa.ref]]
#   seqinr::write.fasta(fa.ref, sub("^>", "", attr(fa.ref, "Annot")),
#                       paste0(tempt.dir, "/mpw.ref.fa"))
#   formatted.header <- sub("^>", "", sapply(fa, attr, "Annot"))  %>% gsub("[^0-9a-zA-Z]", "_",.)
#   seqinr::write.fasta(fa, formatted.header,
#                       paste0(tempt.dir, "/mpw.strain.fa"))

#   align.cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/python3 ",
#                       "/bar/cfan/viralBrowser/release_github_4_22/publicAlignment.py",
#                       " --script_dir ", "/bar/cfan/viralBrowser/release_github_4_22/",
#                       " --ref_fa ", tempt.dir, "/mpw.ref.fa",
#                       " --ref_name ", chr.name,
#                       " --shift ", shift,
#                       " --strain_fa ", tempt.dir, "/mpw.strain.fa",
#                       " --tempt_dir ", tempt.dir,
#                       " --SNV_dir ", SNV.dir, 
#                       # " --email changxu.fan@gmail.com",
#                       " --aligner ", aligner)
#   print(align.cmd)
#   if (dry == F)
#     system(align.cmd)
#   if (out.json==T) {
#     # snvs <- system(paste0("ls ", SNV.dir, "/*bed.gz"), intern = T) %>% basename()
#     jsongen <- lapply(formatted.header, function(x) {
#       df <- data.frame(name = x, url = paste0(x, ".bed.gz"), type = "pairwise")
#     }) %>% Reduce(rbind,.)
#     jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify() %>% write(paste0(SNV.dir, "/snv.json"))
#   }
#   names(formatted.header) <- NULL
#   return(formatted.header)
# }

# readGTF.fanc <- function(gtf.file, get.fields=NULL) {
#   gtf <- read.table(gtf, sep = "\t", quote = "", header = F, as.is = T)
#   colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#   if (!is.null(get.fields))
#     gtf <- gtf.get.fields(gtf, fields = get.fields)
#   return(gtf)
# }

# gtf.get.fields <- function(gtf, fields) {
#   for (field in fields) {
#     gtf[, field] <- str_extract(string = gtf$attribute, pattern = paste0(field,"[^;]+")) %>% gsub("\"", "",.) %>% sub(paste0(field, " "), "",. )
#   }
#   return(gtf)
# }

# # t <- gtf.get.fields(gtf.df, "gene_name")

# stitch.fa <- function(bed, genome, out.file) {
#   # bed has to be written in a way that columns starting from the 4th one represents "group"
#   # sequences are grouped based on all columns (starting from 4th) cat together ("|" separated)
#   tempt.fa <- bed.align(bed = bed, aligned.fa = out.file, genome = genome, add.coordinate = F)

#   fas <- seqinr::read.fasta(tempt.fa)
#   fa.df <- data.frame(header = lapply(fas, attr, "name") %>% unlist() %>% `names<-`(NULL),
#                       seq = lapply(fas, function(x) paste0(x[1:length(x)], collapse = "") ) %>% unlist())
#   fa.df <- fa.df %>% group_by(header) %>% summarise(seq = paste0(seq, collapse = ""))
#   seqinr::write.fasta(sequences = fa.df$seq, names = fa.df$header, file.out = out.file)
#   return(list(fa.df = fa.df, fa.file = fa.file))

# }

# # stitch.fa("~/test/random/test.bed", genome = "mm10", out.file = "~/test/random/test_miao.fa")

# bowtie2.genome.gen <- function(map.bed, genome, bowtie2.build = "/opt/apps/bowtie2/2.3.4.1/bowtie2-build",
#                                out.dir, root.name, thread = 8, run = T) {
#   # map.bed: bed file. only first 3 columns used.currently only supports mono-region.
#   # get fasta file:
#   system(paste0("mkdir -p ", out.dir))
#   map.bed <- map.bed[1, 1:3]
#   map.bed$V4 <- map.bed[, 1]
#   fa <- bed.align(bed = map.bed, aligned.fa = paste0(out.dir, "/", root.name), genome = genome, mask = F, add.coordinate = F)
#   # now make a mini genome:
#   cmd <- paste0("cd ",out.dir," && ", bowtie2.build, " --thread ", thread,
#                 " ", basename(fa), " ", root.name)
#   print(cmd)
#   if (run == T)
#     system(cmd)

#   return(list(bowtie.genome = paste0(out.dir, "/", root.name), map.bed = map.bed, genome = genome))
# }

# mini.align <- function(fastq.vec, R1.pattern="1.fastq.gz", R2.pattern = "2.fastq.gz", k = NULL, a = F,
#                        bowtie.genome, map.bed, out.bam, bowtie.log.file=NULL, thread = 16, run = T,
#                        bowtie2 = "/opt/apps/bowtie2/2.3.4.1/bowtie2", no.unal = F) {
#   system(paste0("mkdir -p ", dirname(out.bam)))

#   R1 <- fastq.vec[grepl(R1.pattern, fastq.vec)][1]
#   R2 <- fastq.vec[grepl(R2.pattern, fastq.vec)][1]
#   # cmd <- paste0("/bar/cfan/scripts/dna-seq/bowtie2_4_22_20.sh ", " -p ", thread,
#   #               " -x ", bowtie.genome, " -i ", R1, " -I ", R2, " -o ", out.bam, " -k ", 50,
#   #               " -s ", map.bed[1,2])
#   cmd <- paste0(bowtie2, " -X2000 --reorder --very-sensitive --xeq --seed 42 ",
#                 " -p ", thread, " -x ", bowtie.genome, " -1 ", R1, " -2 ", R2)

#   # if (!is.null(bowtie.log.file))
#   #   cmd <- paste0(cmd, " --met-file ", bowtie.log.file)

#   if (!is.null(k))
#     cmd <- paste0(cmd, " -k ", k)
#   if (a == T)
#     cmd <- paste0(cmd, " -a ")
#   if (no.unal == T)
#     cmd <- paste0(cmd, " --no-unal")

#   shift <- map.bed[1,2]
#   cmd <- paste0(cmd , " | ", "awk -F \"\\t\" 'BEGIN {OFS = \"\\t\"} {$4 = $4+",shift,"; print $0}' | ")
#   cmd <- paste0(cmd, " /bar/cfan/anaconda2/envs/jupyter/bin/samtools sort -O bam -m 4G -@ ",thread," -  > ", out.bam)

#   utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
#   if (!file.exists(out.bam))
#     stop(paste0(out.bam, " was not successfully generated"))

#   cmd <- paste0("samtools index ", out.bam)
#   utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run, stdout.file = bowtie.log.file)
#   if (!file.exists(out.bam %>% paste0(".bai")))
#     stop(paste0(out.bam, ".bai was not successfully generated"))

#   return(out.bam)
# }

# trimmed.fastq.gen <- function(fastq.vec, thread, run=T) {
#   cmd <- paste0("python3 ", "~/software/atac/encode_pipeline/src/encode_task_trim_adapter.py ", paste0(fastq.vec, collapse = " "),
#                 " --paired-end ", " --auto-detect-adapter ", " --cutadapt-param ' -e 0.1 -m 5' ", " --nth ", thread)
#   utilsFanc::cmd.exec.fanc(cmd = cmd, run = run, intern = F)
# }


align.by.gtf <- function(abao, gtf.gr, exon_numbers = NULL,
                         key.gene, key.regionID.in.abao, region.type = "CDS", buffer.left = 0, buffer.right = 0,
                         regionIDs.exclude = NULL,
                         out.dir) {
  # specify buffer as if you are expanding the region. left = 20, right = 20 will expand 40, 20 on each side
  # negative integers will shrink the boundaries
  gr <- gtf.gr %>% plyranges::filter(gene_name == key.gene, type == region.type)
  if (!is.null(exon_numbers)) {
    gr <- gr %>% plyranges::filter(exon_number %in% as.character(exon_numbers))
  }
  if (any(duplicated(gr$exon_number))) {
    stop("any(duplicated(gr$exon_number))")
  }
  system(paste0("mkdir -p ", out.dir, "/smart_cut/"))
  abl <- lapply(1:length(gr), function(i) {
    region <- gr[i]
    int.name <- paste0(region.type, "_", region$exon_number)
    if (buffer.left > 0 && buffer.left < 1) {
      buffer.left <- floor(buffer.left * width(region))
    }
    if (buffer.right > 0 && buffer.right < 1) {
      buffer.right <- floor(buffer.right * width(region))
    }
    start(region) <- start(region) - buffer.left
    end(region) <- end(region) + buffer.right
    
    smart.cut <- region %>% `names<-`(NULL) %>% as.data.frame() %>% 
      dplyr::select(seqnames, start, end) %>% dplyr::rename(chr = seqnames) %>% 
      dplyr::mutate(int.name = int.name, regionID = key.regionID.in.abao, 
             buffer.left = 0, buffer.right = 0)
    smart.cut.tsv <- paste0(out.dir, "/smart_cut/", int.name, ".tsv")
    write.table(smart.cut, smart.cut.tsv, quote = F, sep = "\t", row.names = F, col.names = T)
    abao.sub <- aba.subset(abao = abao, smart.cut.df = smart.cut.tsv,
                           regionIDs.exclude = regionIDs.exclude,
                           new.work.dir = paste0(out.dir, "/", int.name, "/"), 
                           threads = 1)
    return(abao.sub)
  })
  return(abl)
}



percent.identity.matrix <- function(aligned.fa = NULL, start = 1, end = NULL, 
                                    seq.vec, ident.mat = NULL, PID = "PID1", 
                                    subset.rows = NULL, subset.cols = NULL,
                                    color.range = c(0.5, 1), label.hm = T,
                                    width = 2100, height = 2000,
                                    out.dir = NULL, root.name, 
                                    do.save.rds = T, do.plot = T,
                                    # save.rds = NULL, plot.out = NULL,
                                    write.pairwise.aln = F,
                                    threads = 1, ...) {
  if (do.save.rds || do.plot) {
    if (is.null(out.dir)) {
      out.dir <- dirname(aligned.fa)
    }
    if (is.null(root.name)) {
      root.name <- tools::file_path_sans_ext(basename(aligned.fa))
    }
  }
  
  if (is.null(ident.mat)) {
    if (!is.null(aligned.fa)) {
      aln <- readDNAMultipleAlignment(aligned.fa)
      browser()
      seq.vec <- as(aln,"DNAStringSet") %>% as.character()
    }
    seq.vec <- seq.vec %>% unlist()
    if (write.pairwise.aln == T) {
      aln.dir <- paste0(out.dir, "/", root.name, "/")
    } else {
      aln.dir <- NULL
    }
    ident.mat <- outer(seq.vec, seq.vec, function(x, y) 
      pairwise.identity(seq1s = x, seq2s = y, PID = PID,
                        threads = threads, aln.out.dir = aln.dir, ...))

    if (do.save.rds) {
      rds <- paste0(out.dir, "/", root.name, "_ident_mat.Rds")
      dir.create(dirname(rds), recursive = T, showWarnings = F)
      saveRDS(ident.mat, rds)
    }
    
  }
  
  rows <- rownames(ident.mat)
  if (!is.null(subset.rows)) {
    rows <- rows %>% .[grepl(subset.rows,.)]
  } 
  
  cols <- colnames(ident.mat)
  if (!is.null(subset.cols)) {
    if (subset.cols[1] == "rows") {
      cols <- rows
    } else {
      cols <- cols %>% .[grepl(subset.cols,.)]
    }
  }
  ident.mat <- ident.mat[rows, cols]
  
  if (do.plot) {
    plot.out <- paste0(out.dir, "/", root.name,"_ident_mat.pdf")
    plot <- v4c::contact.map.viz(mat = ident.mat, color.range = color.range, 
                                 add.text = label.hm, digits = 3, plot.out = plot.out,
                                 width = width, height = height)
  }
  
  return(ident.mat)
}

# pairwise.identity.fast <- function(seq1s, seq2s, threads = 1) {
#   # assumming they are coming out of a pairwise alignment
#   if (length(seq1s) != length(seq2s)) {
#     stop("length(seq1s) != length(seq2s)")
#   }
#   if (!is.character(seq1s) || !is.character(seq2s)) {
#     stop("!is.character(seq1s) || !is.character(seq2s)")
#   }
#   if (!identical(nchar(seq1s), nchar(seq2s))) {
#     stop("!identical(nchar(seq1s), nchar(seq2s))")
#   }
#   
#   res <- utilsFanc::safelapply(seq_along(seq1s), function(i) {
#     seq1 <- seq1s[[i]]
#     seq2 <- seq2s[[i]] 
#     # couldn't figure out how to get nMatch...
#     
#     return(pct)
#   }) %>% unlist()
#   return(res)
# }
pairwise.identity <- function(seq1s, seq2s, PID = "PID1", threads = 1, 
                              aln.out.dir = NULL, root.name = NULL, ...) {
  if (length(seq1s) != length(seq2s)) {
    stop("length(seq1s) != length(seq2s)")
  }
  if (is.null(names(seq1s))) {
    names(seq1s) <- paste0("seq1_", 1:length(seq1s))
  }
  if (is.null(names(seq2s))) {
    names(seq2s) <- paste0("seq2_", 1:length(seq2s))
  }
  res <- utilsFanc::safelapply(seq_along(seq1s), function(i) {
    seq1 <- seq1s[[i]] %>% DNAstring.rm.gap()
    seq2 <- seq2s[[i]] %>% DNAstring.rm.gap()
    
    paln <- pairwiseAlignment(pattern = seq1, subject = seq2, ...)
    if (!is.null(aln.out.dir)) {
      dir.create(aln.out.dir, showWarnings = F, recursive = T)
      if (is.null(root.name)) {
        root.name <- basename(aln.out.dir)
      }
      name1 <- names(seq1s)[i]
      name2 <- names(seq2s)[i]
      out.fa <- paste0(aln.out.dir, "/", root.name, "_", name1, "-", name2, ".txt")
      writePairwiseAlignments(paln, out.fa)
    }
    pct <- pid(paln, PID) 
    pct <- round(pct/100, digits = 4)
    return(pct)
  }, threads = threads) %>% unlist()
  return(res)
}

DNAstring.rm.gap <- function(x) {
  if (! is.character(x))
    x <- x %>% as.character()
  x %>% gsub("-+", "", .) %>% DNAString() %>% return()
}

pairwise.ident.slide <- function(aligned.fa, start = 1, end, step.size = 20, window.size = 20) {
  start.vec <- seq(from = start, to = end-window.size - 1, by = step.size)
  end.vec <- start.vec + step.size
  browser()
}

pairwise.dist.mat <- function(input.list, fg.grep = "fg", comps, comp.is.file = T, slide.fg = F, model = "TN93",
                              win.size = 250, slide.size = 250, 
                              last.window.merge.threshold = 100,
                              out.file = NULL) {
  
  # comps: list(acti = list(clade1 = c("Ly49d", "Ly49m"), clade2 = c("Ly49h", "Ly49n")),
  #             inhi = list(clade1 = c("Ly49a", "Ly49g"), clade2 = c("Ly49c", "Ly49i")))
  #comps serves as 2 filters: 1. only genes included in comps will be used. 2. only inter clade 
  #distances are kept unless there is only 1 clade. clade names are not important.
  # the final estimate of will be the mean of **inter** clade distances
  # if (any(!grepl("^fg|^bg", names(input.list)))) {
  #   stop("input.list name error")
  # }
  input.list <- lapply(input.list, function(x) {
    return(Biostrings::readDNAMultipleAlignment(x) %>% as.matrix())
  })
  input.genes <- Reduce(intersect, lapply(input.list, rownames))
  stats <- lapply(names(comps), function(comp.name) {
    comp <- comps[[comp.name]]
    stats <- lapply(names(input.list), function(input.name) {
      seq.mat <- input.list[[input.name]][input.genes, ]
      bFg <- grepl(paste0(fg.grep, collapse = "|"), input.name)
      
      #### start sliding window
      if (bFg == F || slide.fg == T) {
        df <- data.frame(chr = "chr1", start = 1, end = ncol(seq.mat))
        gr <- GenomicRanges::slidingWindows(x = makeGRangesFromDataFrame(df),
                                            width = win.size, step = slide.size)
        slide.df <- gr %>% as.data.frame() %>% dplyr::select(start, end)
        # note: somehow slidingwindows return a grange list object.
        if (width(gr)[[1]][length(gr[[1]])] < last.window.merge.threshold) {
          slide.df[nrow(slide.df) - 1, "end"] <- slide.df[nrow(slide.df), "end"]
          slide.df <- slide.df[-nrow(slide.df), ]
        }
      } else {
        slide.df <- data.frame(start = 1, end = ncol(seq.mat))
      }
      
      slide.df$id <- paste0(slide.df$start, "_", slide.df$end)
      
      stats <- lapply(1:nrow(slide.df), function(i) {
        start <- slide.df[i, "start"]
        end <- slide.df[i, "end"]
        stats <- pairwise.dist.mat.core(seq.mat = seq.mat[,start:end], comp = comp,
                                        comp.is.file = comp.is.file,
                                        model = model)$stats
        return(stats)
      }) %>% Reduce(rbind, .)
      stats <- utilsFanc::add.column.fanc(df1 = stats,
                                          df2 = data.frame(comp = comp.name, input = input.name, 
                                                           window = slide.df$id),
                                          pos = 1)
      return(stats)
    }) %>% Reduce(rbind, .)
    return(stats)
  }) %>% Reduce(rbind, .)
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), recursive = T, showWarnings = F)
    write.table(stats, out.file, sep = "\t", row.names = F, col.names = T, quote = F)
  }
  return(stats)
}

pairwise.dist.mat.core <- function(seq.mat, model = "TN93", comp = NULL, comp.is.file = F,
                                   use.pct.ident = F, pairwise.deletion = T,
                                   plot.out = NULL, height = 500, width = 600, 
                                   color.range = c(0, 0.6)) {
  if (file.exists(seq.mat[1])) {
    seq.mat <- Biostrings::readDNAMultipleAlignment(filepath = seq.mat) %>% 
      Biostrings::as.matrix()
  }
  if (comp.is.file) {
    comp <- lapply(comp, readLines)
  }
  if (is.null(comp)) {
    comp <- list(ungroupped = rownames(seq.mat))
  }

  comp.genes <- comp %>% unlist()
  if (any(duplicated(comp.genes))) {
    stop("any(duplicated(comp.genes))")
  }
  genes.use <- comp.genes %>% .[. %in% rownames(seq.mat)]
  if (length(genes.use) < 2) {
    stop("length(genes.use) < 2")
  }
  seq.mat <- seq.mat[genes.use, ]
  dist.mat <- ape::dist.dna(x = ape::as.DNAbin(seq.mat), as.matrix = T, model = model,
                            pairwise.deletion = pairwise.deletion)
  if (use.pct.ident) {
    dist.mat <- 1 - dist.mat
    color.range <- sort(1-color.range)
  }
  name.mat <- outer(rownames(dist.mat), rownames(dist.mat), paste, sep = "-")
  groups <- names(comp)
  if (length(groups) > 1) {
    mapping.vec <- lapply(groups, function(group) {
      x <- rep(group, length(comp[[group]]))
      names(x) <- comp[[group]]
      return(x)
    }) %>% `names<-`(NULL) %>% Reduce(c, .)
    g.vec <- mapping.vec[rownames(dist.mat)]
    paste.mat <- outer(g.vec, g.vec, FUN = paste0)
    mask.values <- paste0(groups, groups)
    bCross <- !paste.mat %in% mask.values
    dist.vec <- dist.mat[bCross]
    names(dist.vec) <- name.mat[bCross]
  } else if (length(groups) == 1) {
    dist.vec <- dist.mat %>% as.vector()
    names(dist.vec) <- dist.mat %>% as.vector()
  } else {
    stop("length(groups) < 1")
  }
  stats <- data.frame(mean = mean(dist.vec), sd = stats::sd(dist.vec))
  if (!is.null(plot.out)) {
    v4c::contact.map.viz(mat = dist.mat,add.text = T, digits = 3,
                         height = height, width = width,
                         color.range = color.range, 
                         plot.out = plot.out)
  }
  return(list(dist.mat = dist.mat, dist.vec = dist.vec, 
              stats = stats))
}

pairwise.dist.plot <- function(stats, meta.df = NULL, linetype, fg.grep = "fg",
                               add.boundary = T, width = 8, height = 4,
                               pt.size = 0.1, line.size = 1,
                               plot.out = NULL) {
  if (is.character(stats)) {
    stats <- read.table(stats, header = T)
  }
  if (!is.null(meta.df)) {
    stats <- left_join(stats,meta.df)
  }
  stats$region <- paste0(stats$input, ":", stats$window)
  stats$region <- factor(stats$region, levels = unique(stats$region))
  stats$type <- "bg"
  stats$type[grepl(paste0(fg.grep, collapse = "|"), stats$input)] <- "fg"
  
  p <- ggplot(stats, aes(x = region, y = mean)) + 
    geom_point(aes(shape = type), size = pt.size)
  if (is.null(meta.df)) {
    p <- p + geom_line(aes(group = comp, color = comp), size = line.size)
  } else {
    p <- p + geom_line(aes_string(group = "comp", color = "comp", linetype = linetype), size = line.size)
  }
  p <- p+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (add.boundary == T) {
    boundaries <- stats %>% group_by(input) %>% summarise(start = window[1]) %>% 
      dplyr::mutate(bound = paste0(input, ":", start)) %>% pull(bound)
    p <- p + geom_vline(xintercept = boundaries)
  }
  if (!is.null(plot.out))
    scFanc::wrap.plots.fanc(list(p), sub.height = height, sub.width = width, 
                            plot.out = plot.out)
  invisible(p)
}

fa.cons.shrink <- function(fa, out.fa = NULL, max.gap.frac = 0.3,
                           regionIDs.include = NULL, regionIDs.exclude = NULL,
                           # parameters copied from aba.add.consensus
                           ambiguityMap = "N", threshold = 0.5, regions.use = NULL,
                           consensus.name = NULL, consensus.doc = NULL) {
  seq.list <- seqinr::read.fasta(fa)

  if (is.null(consensus.name)) {
    if (identical(ambiguityMap, Biostrings::IUPAC_CODE_MAP))
      map.name <- "IUPAC"
    else
      map.name <- ambiguityMap
    consensus.name <- paste0("cons_", map.name, "_", threshold)
  }
  if (is.null(consensus.doc)) {
    consensus.doc <- consensus.name
  }
  
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
  cons.score <- msa::msaConservationScore(cons.mat, sub.mat, gapVsGap=0)
  cons.string <- Biostrings::consensusString(strings, ambiguityMap = ambiguityMap, 
                                             threshold = threshold)
  shrink.map <- data.frame(seq = cons.string %>% strsplit(split = "") %>% unlist(),
                           pos.full = 1:nchar(cons.string))
  shrink.map <- shrink.map %>% filter(seq != "-")
  shrink.map$pos.shrink <- 1:nrow(shrink.map)
  seq.list.shrink <- lapply(seq.list, function(x) {
    return(x[shrink.map$pos.full])
  })
  cons.list <- list(in.fa = fa,
                    seq = cons.string, 
                    seq.shrink = cons.string %>% gsub("-+", "", .),
                    cons.score = cons.score,
                    cons.score.shrink = cons.score[names(cons.score) != "-"],
                    shrink.map = shrink.map,
                    name = consensus.name,
                    doc = consensus.doc)
  if (is.null(out.fa)) {
    out.fa <- paste0(dirname(fa), "/", tools::file_path_sans_ext(basename(fa)), "_shrink.fa")
  }
  aba.write.consensus.with.aln.2(aln.list = seq.list.shrink, max.gap.frac = max.gap.frac,
                                 out.file = out.fa, regionIDs.include = regionIDs.include,
                                 regionIDs.exclude = regionIDs.exclude)
  out.Rds <- paste0(tools::file_path_sans_ext(out.fa), ".Rds")
  saveRDS(cons.list, out.Rds)
  return(out.fa)
}


pairwise.dist.mat.core <- function(seq.mat, model = "TN93", comp = NULL, comp.is.file = F,
                                   use.pct.ident = F, pairwise.deletion = T,
                                   plot.out = NULL, height = 500, width = 600, 
                                   color.range = c(0, 0.6)) {
  if (file.exists(seq.mat[1])) {
    seq.mat <- Biostrings::readDNAMultipleAlignment(filepath = seq.mat) %>% 
      Biostrings::as.matrix()
  }
  if (comp.is.file) {
    comp <- lapply(comp, readLines)
  }
  if (is.null(comp)) {
    comp <- list(ungroupped = rownames(seq.mat))
  }

  comp.genes <- comp %>% unlist()
  if (any(duplicated(comp.genes))) {
    stop("any(duplicated(comp.genes))")
  }
  genes.use <- comp.genes %>% .[. %in% rownames(seq.mat)]
  if (length(genes.use) < 2) {
    stop("length(genes.use) < 2")
  }
  seq.mat <- seq.mat[genes.use, ]
  dist.mat <- ape::dist.dna(x = ape::as.DNAbin(seq.mat), as.matrix = T, model = model,
                            pairwise.deletion = pairwise.deletion)
  if (use.pct.ident) {
    dist.mat <- 1 - dist.mat
    color.range <- sort(1-color.range)
  }
  name.mat <- outer(rownames(dist.mat), rownames(dist.mat), paste, sep = "-")
  groups <- names(comp)
  if (length(groups) > 1) {
    mapping.vec <- lapply(groups, function(group) {
      x <- rep(group, length(comp[[group]]))
      names(x) <- comp[[group]]
      return(x)
    }) %>% `names<-`(NULL) %>% Reduce(c, .)
    g.vec <- mapping.vec[rownames(dist.mat)]
    paste.mat <- outer(g.vec, g.vec, FUN = paste0)
    mask.values <- paste0(groups, groups)
    bCross <- !paste.mat %in% mask.values
    dist.vec <- dist.mat[bCross]
    names(dist.vec) <- name.mat[bCross]
  } else if (length(groups) == 1) {
    dist.vec <- dist.mat %>% as.vector()
    names(dist.vec) <- dist.mat %>% as.vector()
  } else {
    stop("length(groups) < 1")
  }
  stats <- data.frame(mean = mean(dist.vec), sd = stats::sd(dist.vec))
  if (!is.null(plot.out)) {
    v4c::contact.map.viz(mat = dist.mat,add.text = T, digits = 3,
                         height = height, width = width,
                         color.range = color.range, 
                         plot.out = plot.out)
  }
  return(list(dist.mat = dist.mat, dist.vec = dist.vec, 
              stats = stats))
}

pairwise.dist.plot <- function(stats, meta.df = NULL, linetype, fg.grep = "fg",
                               add.boundary = T, width = 8, height = 4,
                               pt.size = 0.1, line.size = 1,
                               plot.out = NULL) {
  if (is.character(stats)) {
    stats <- read.table(stats, header = T)
  }
  if (!is.null(meta.df)) {
    stats <- left_join(stats,meta.df)
  }
  stats$region <- paste0(stats$input, ":", stats$window)
  stats$region <- factor(stats$region, levels = unique(stats$region))
  stats$type <- "bg"
  stats$type[grepl(paste0(fg.grep, collapse = "|"), stats$input)] <- "fg"
  
  p <- ggplot(stats, aes(x = region, y = mean)) + 
    geom_point(aes(shape = type), size = pt.size)
  if (is.null(meta.df)) {
    p <- p + geom_line(aes(group = comp, color = comp), size = line.size)
  } else {
    p <- p + geom_line(aes_string(group = "comp", color = "comp", linetype = linetype), size = line.size)
  }
  p <- p+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (add.boundary == T) {
    boundaries <- stats %>% group_by(input) %>% summarise(start = window[1]) %>% 
      dplyr::mutate(bound = paste0(input, ":", start)) %>% pull(bound)
    p <- p + geom_vline(xintercept = boundaries)
  }
  if (!is.null(plot.out))
    scFanc::wrap.plots.fanc(list(p), sub.height = height, sub.width = width, 
                            plot.out = plot.out)
  invisible(p)
}

fa.cons.shrink <- function(fa, out.fa = NULL, max.gap.frac = 0.3,
                           regionIDs.include = NULL, regionIDs.exclude = NULL,
                           # parameters copied from aba.add.consensus
                           ambiguityMap = "N", threshold = 0.5, regions.use = NULL,
                           consensus.name = NULL, consensus.doc = NULL) {
  seq.list <- seqinr::read.fasta(fa)

  if (is.null(consensus.name)) {
    if (identical(ambiguityMap, Biostrings::IUPAC_CODE_MAP))
      map.name <- "IUPAC"
    else
      map.name <- ambiguityMap
    consensus.name <- paste0("cons_", map.name, "_", threshold)
  }
  if (is.null(consensus.doc)) {
    consensus.doc <- consensus.name
  }
  
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
  cons.score <- msa::msaConservationScore(cons.mat, sub.mat, gapVsGap=0)
  cons.string <- Biostrings::consensusString(strings, ambiguityMap = ambiguityMap, 
                                             threshold = threshold)
  shrink.map <- data.frame(seq = cons.string %>% strsplit(split = "") %>% unlist(),
                           pos.full = 1:nchar(cons.string))
  shrink.map <- shrink.map %>% filter(seq != "-")
  shrink.map$pos.shrink <- 1:nrow(shrink.map)
  seq.list.shrink <- lapply(seq.list, function(x) {
    return(x[shrink.map$pos.full])
  })
  cons.list <- list(in.fa = fa,
                    seq = cons.string, 
                    seq.shrink = cons.string %>% gsub("-+", "", .),
                    cons.score = cons.score,
                    cons.score.shrink = cons.score[names(cons.score) != "-"],
                    shrink.map = shrink.map,
                    name = consensus.name,
                    doc = consensus.doc)
  if (is.null(out.fa)) {
    out.fa <- paste0(dirname(fa), "/", tools::file_path_sans_ext(basename(fa)), "_shrink.fa")
  }
  aba.write.consensus.with.aln.2(aln.list = seq.list.shrink, max.gap.frac = max.gap.frac,
                                 out.file = out.fa, regionIDs.include = regionIDs.include,
                                 regionIDs.exclude = regionIDs.exclude)
  out.Rds <- paste0(tools::file_path_sans_ext(out.fa), ".Rds")
  saveRDS(cons.list, out.Rds)
  return(out.fa)
}



pairwise.dist.mat.ez <- function(seq.mat, model, genes = NULL, regex.genes = F,
                                   use.pct.ident = F, pairwise.deletion = T,
                                   plot.out = NULL, plot.dir = NULL, suffix = NULL, skip.plot = F,
                                 # by default: modify the input file name.
                                 height = 500, width = 600, plot.add.text = T,
                                 plot.cluster.seqs = F, plot.cluster.method = "euclidean",
                                   color.range = c(0, 0.6)) {
  if (file.exists(seq.mat[1])) {
    if (is.null(plot.out)) {
      if (is.null(plot.dir)) {
        plot.dir <- paste0(tools::file_path_sans_ext(seq.mat), "_dist/")
      }
      if (is.null(suffix)) {
        suffix <- Sys.time() %>% as.character() %>% gsub(" +", "_", .)
      }
      root <- basename(tools::file_path_sans_ext(seq.mat))
      plot.out <- paste0(plot.dir, "/", root, "_", suffix, ".png")
    }
    seq.mat <- Biostrings::readDNAMultipleAlignment(filepath = seq.mat) %>% 
      Biostrings::as.matrix()
  }
  if (is.null(rownames(seq.mat))) {
    stop("is.null(rownames(seq.mat))")
  }
  if (is.null(genes)) {
    genes <- rownames(seq.mat)
  }
  if (regex.genes) {
    genes <- rownames(seq.mat) %>% .[grepl(paste0(genes, collapse = "|"), .)]
  }
  genes <- genes %>% .[.%in% rownames(seq.mat)]
  if (length(genes) < 1) {
    stop("length(genes) < 1")
  }
  seq.mat <- seq.mat[genes, ]
  dist.mat <- ape::dist.dna(x = ape::as.DNAbin(seq.mat), as.matrix = T, model = model,
                            pairwise.deletion = pairwise.deletion)
  if (use.pct.ident) {
    dist.mat <- 1 - dist.mat
    color.range <- sort(1-color.range)
  }
  
  if (!is.null(plot.out) && !skip.plot) {
    v4c::contact.map.viz(mat = dist.mat,add.text = plot.add.text, digits = 3,
                         height = height, width = width, 
                         cluster_rows = plot.cluster.seqs,
                         cluster_columns = plot.cluster.seqs,
                         clustering_distance_rows = plot.cluster.method,
                         clustering_distance_columns = plot.cluster.method,
                         color.range = color.range, 
                         plot.out = plot.out)
  }
  return(dist.mat)
}

# pairwise.dist.random.simu <- function(
#   n.pairs = 1000, seed = 42,seq.length = 1000,
#   genome = "mm10", seq,
#   threads = 12,
#   out.dir, root = NULL) {
#   if (is.null(root)) {
#     root <- basename(out.dir)
#   }
#   set.seed(seed = seed)
#   
#   utilsFanc::safelapply(1:n.pairs, function(i) {
#     
#   })
#   
# }

slide.pid <- function(abao, ref.id, seqs.use.regex = NULL, seqs.rev.regex = NULL,
                      win.size = 100, conservation.pids = 70, conservation.length = 100,
                      align.method = "mafft", pid.method = "vista",
                      regionID.color.map = NULL, pid.color.map = NULL,
                      ymin = 50, ymax = 85, remove.below.ymin = F,
                      npar.query = 1, rgbPeak.genome = "mm10",
                      out.dir = NULL, root.name = NULL, ...) {
  abao <- aba.import(abao)
  rownames(abao@ori.df) <- abao@ori.df$regionID
  fa <- paste0(abao@work.dir, "/", abao@meta.data$fa.root.name, "in.fa")
  seq.list <- seqinr::read.fasta(fa, forceDNAtolower = F)
  if (!ref.id %in% names(seq.list)) {
    stop("!ref.id %in% names(seq.list)")
  }
  
  if (!is.null(seqs.rev.regex)) {
    seqs.rev <- names(seq.list) %>% .[grepl(paste0(seqs.rev.regex, collapse = "|"), .)]
    print(paste0("reverse complementing these sequences: ", 
                 paste0(seqs.rev, collapse = ", ")))
    seq.list[seqs.rev] <- lapply(seq.list[seqs.rev], function(x) {
      return(rev(seqinr::comp(x, forceToLower = F)))
    })
  }
  
  if (!is.null(seqs.use.regex)) {
    seqs.query <- names(seq.list) %>% .[grepl(paste0(seqs.use.regex, collapse = "|"), .)]
  } else {
    seqs.query <- names(seq.list)
  }
  seqs.query <- seqs.query %>% .[.!= ref.id]
  if (length(seqs.query) < 1) {
    stop("length(seqs.query) < 1")
  }
  seq.list <- seq.list[c(ref.id, seqs.query)]
  
  if (is.null(out.dir)) {
    # out.dir <- fa %>% tools::file_path_sans_ext() %>% paste0("_slidePid/")
    out.dir <- paste0(abao@work.dir, "/slidePid/ref", ref.id, "/")
  }
  if (is.null(root.name)) {
    # root.name <- paste0(abao@meta.data$fa.root.name, "slide_ref", ref.id)
    root.name <- abao@meta.data$fa.root.name %>% sub("_$", "", .)
  }
  pair.dir <- paste0(out.dir, "/pairwise/")
  # seq.dir <- paste0(out.dir, "/seqs/")
  dir.create(pair.dir, showWarnings = F, recursive = T)
  # dir.create(seq.dir, showWarnings = F, recursive = T)
  
  if (align.method == "mafft") {
    pairs <- utilsFanc::safelapply(seqs.query, function(qname) {
      aln.fa <- paste0(pair.dir, "/", root.name, "_", qname, "_vs_", ref.id,
                       "_", align.method, ".fa")
      trash <- mafft.fanc(in.fa = seq.list[c(ref.id, qname)], aligned.fa = aln.fa, 
                          mafft = MAFFT.DEFAULT, mafft.params = MAFFT.AUTO)
      return(aln.fa)
    }, threads = npar.query) %>% unlist()
    names(pairs) <- seqs.query
  }
  if (tolower(pid.method) == "vista") {
    vista.in <- utilsFanc::safelapply(pairs, fa.2.vista, ref.id = ref.id, threads = npar.query) %>% unlist()
    names(vista.in) <- names(pairs)
    vista.dir <- paste0(out.dir, "/vista/")
    vista.root <- paste0(root.name, "_", ref.id)
    lapply(conservation.pids, function(conservation.pid) {
      vista.run(vistain.vec = vista.in, ref.name = ref.id,
                root.name = vista.root, out.dir = vista.dir, run = T,
                win.size =win.size, conservation.pid = conservation.pid,  
                conservation.length = conservation.length, ...)
      
    })
    # conservation bedgraph's and conservation bed's:
    browser.df <- lapply(names(vista.in), function(q) {
      bdg <- vista.score.2.bdg(score.file = paste0(vista.dir, "/", vista.root, "_", q, "_scores.txt"), 
                        chr = abao@ori.df[ref.id, "chr"], start = abao@ori.df[ref.id, "start"],
                        min.show = ifelse(remove.below.ymin, ymin, 0))
      beds <- sapply(conservation.pids, function(conservation.pid) {
        region.file <- paste0(vista.dir, "/", vista.root, "_", q, "_cons_", conservation.pid, "_", 
                              conservation.length,"_regions.txt")
        bed <- vista.region.2.bed(region.file = region.file, 
                                  chr = abao@ori.df[ref.id, "chr"], start = abao@ori.df[ref.id, "start"],
                                  must.exist = F)
        return(bed)
      })
      
      if (!is.null(pid.color.map)) {
        rownames(pid.color.map) <- as.character(pid.color.map$pid)
        pid.colors <- pid.color.map[as.character(conservation.pids), "color"]
        if (any(is.na(pid.colors))) {
          stop("any(is.na(pid.colors))")
        }
      } else {
        pid.colors <- NULL
      }
      # if (q == "dog.Ly49") {
      #   browser()
      # }
      rgbPeak <- paste0(vista.dir, "/", vista.root, "_", q, "_cons_", 
                        paste0(conservation.pids, collapse = "_"), "_", 
                        conservation.length,"_regions.rgbPeak")
      rgbPeak <- utilsFanc::merge.bed.color(beds = beds, colors = pid.colors, 
                                 out.rgbPeak = rgbPeak, genome = rgbPeak.genome)
      
      # browser.df <- data.frame(regionID = q, name = paste0(q, "_", c("0pid", paste0(1:length(conservation.pids), "cons_", conservation.pids))), 
      #                          url = paste0(utilsFanc::bash2ftp(c(bdg, beds)), ".gz"), 
      #                          track_type = c("bedgraph", rep("bed", length(beds))),
      #                          pid = c(NA, conservation.pids))
      # # 0pid: 0 makes sure that pid track is on top of the cons bed tracks
      browser.df <- data.frame(regionID = q, name = paste0(q, "_", c("pid", paste0("cons_", paste0(conservation.pids, collapse = "_")))), 
                               url = utilsFanc::bash2ftp(c(paste0(bdg, ".gz"), rgbPeak)),
                               track_type = c("bedgraph", "rgbPeak"))
      return(browser.df)
    }) %>% do.call(rbind, .)
    browser.df$metadata <- "novalue"
    browser.df$options <- ""
    browser.df$options[browser.df$track_type == "bedgraph"] <- paste0("yScale:fixed|yMin:", ymin, "|", "yMax:", ymax)
    browser.df$options[browser.df$track_type == "rgbPeak"] <- "height:10"
    browser.dfs <- browser.df %>% split(., f = factor(.$track_type, levels = c("bedgraph", "rgbPeak")))
    
    sub.f.add.color <- function(track.df, join.by) {
      if (join.by == "regionID") {
        map.df <- regionID.color.map
      } else if (join.by == "pid") {
        map.df <- pid.color.map
      }
      
      if (!is.null(map.df)) {
        if (is.character(map.df)) map.df <- read.table(map.df, header = T)
        utilsFanc::check.intersect(c(join.by, "color"), "required columns", 
                                   colnames(map.df), "colnames of map.df")
        map.df <- map.df[, c(join.by, "color")]
        track.df <- dplyr::left_join(track.df, map.df, by = join.by)
        bNA <- is.na(track.df$color)
        track.df$color[!bNA] <- paste0("color:", track.df$color[!bNA])
        track.df$color[bNA] <- "novalue"
        track.df <- track.df %>% 
          dplyr::mutate(options = paste0(options, "|", color), color = NULL)
        
      }
      return(track.df)
    }
    browser.dfs$bedgraph <- sub.f.add.color(browser.dfs$bedgraph, "regionID")
    # browser.dfs$bed <- sub.f.add.color(browser.dfs$bed, "pid")
    browser.df <- do.call(rbind, browser.dfs)
    browser.df$options <- browser.df$options %>% gsub("\\|+" , "|", .) %>% 
      sub("^\\|", "", .) %>% sub("\\|$", "", .)
    browser.df$options[browser.df$options == ""] <- "novalue"
    browser.df$regionID <- factor(browser.df$regionID, levels = seqs.query)
    browser.df$track_type <- factor(browser.df$track_type, levels = c("bedgraph", "rgbPeak"))
    browser.df <- browser.df %>% dplyr::arrange(regionID, track_type)
    
    json.out <- paste0(out.dir, "/", root.name, "_", ref.id, ".json")
    json <- utilsFanc::jsongen(
      df = browser.df, 
      outfile = json.out)
    cat(utilsFanc::bash2ftp(json.out))
    cat("\n")
  }
  return(json.out)
}


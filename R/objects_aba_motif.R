aba.label.fimo <- function(abao, fimo.dir, p.cutoff = 0.0001, is.hocomoco = F, 
                           is.shrink, 
                           motifs = NULL, arche.exclude = c("GC-tract"), 
                           motif.tf.map = "~/motifs/archetype/motif_tf.tsv",
                           gex.df = NK.GEX, gex.cutoff = 30, 
                           collapse = F,
                           rank.plot = F, to.bed = F, out.dir, root.name = NULL, # out.dir and root.name are used to plot gex rank, if to.bed = T.
                           split.by.int.name = F, split.by = "genome") {
  # fimo can't scan through gaps. So shrunk cons is usually used. 
  # cons.name no longer needed: it's written in fimo's report
  # warning("note: the strand info on the generated bed files reflect the strand of the sequence in ori.df, not the motif strand")
  # warning canceled: I disabled strand info in bed. the strand info in aln.df is directly from fimo
  warning("the flaw of this function: in returned aln.df, all starts and ends are of the coordinates of aln" %>% 
            paste0(" regardless of chr"))
  fimo.df <- read.table(paste0(fimo.dir, "/fimo.tsv"), header = T, sep = "\t")
  fimo.df <- fimo.df %>% filter(p.value <= p.cutoff)
  if (is.hocomoco == T)
    fimo.df$motif_id <- fimo.df$motif_id %>% sub("_.+$", "", .)
  else
    stop("currently only developed for hocomoco motifs")
  if (!is.null(motifs))
    fimo.df <- fimo.df %>% filter(motif_id %in% motifs)
  
  aln.df <- fimo.df %>% mutate(chr = sequence_name, int.name = motif_id, end = stop) %>% 
    select(chr, start, end, int.name, strand)
  
  aln.df$gene <- aln.df$int.name %>% stringr::str_to_title()
  aln.df <- aln.df %>% left_join(gex.df)
  aln.df <- aln.df %>% na.omit()
  
  aln.df <- aln.df %>% split(., f = .$chr) %>% 
    lapply(function(df) {
      chr <- df$chr[1]
      if (chr %in% names(abao@meta.data$consensus)) {
        if (is.shrink == T) {
          df$start <- aba.shrink.2.full(abao = abao, pos = df$start, cons.name = chr)
          df$end <- aba.shrink.2.full(abao = abao, pos = df$end, cons.name = chr)
        }
        return(df)
      } else if (chr %in% abao@ori.df$regionID) {
        df$start <- aba.map.core(abao = abao, regionID = chr, pos = df$start)
        df$end <- aba.map.core(abao = abao, regionID = chr, pos = df$end)
        return(df)
      } else {
        return(NULL)
      }
    }) %>% Reduce(rbind, .)
  
  p <- scFanc::rank.plot(df = gex.df, vars = "mean", transformation = function(x) log2(x + 1), label.size = 2,
                         rank.method = "random", label.var = "gene", add.h.line = gex.cutoff,
                         labels = aln.df$gene, outfile = paste0(out.dir, "/", root.name, "_gex_rank.png"))
  
  aln.df <- aln.df %>% filter(mean > gex.cutoff)
  if (is.character(motif.tf.map))
    motif.tf.map <- read.table(motif.tf.map, header = T)
  
  motif.tf.map$gene <- motif.tf.map$gene %>% stringr::str_to_title()
  motif.tf.map <- motif.tf.map %>% .[!duplicated(.$gene),]
  motif.tf.map <- utilsFanc::change.name.fanc(df = motif.tf.map, cols.from = "Name", cols.to = "arche.name")
  
  aln.df <- left_join(aln.df, motif.tf.map) %>% na.omit()
  aln.df <- aln.df %>% arrange(chr, start)
  aln.df$arche.name <- aln.df$arche.name %>% sub("_.+$", "", .)
  if (collapse == T) {
    aln.df <- aln.df %>% split(., f= .$arche.name) %>% 
      lapply(function(df) {
        gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = T)
        gr <- plyranges::reduce_ranges(gr, arche.name = unique(arche.name))
        out <- gr %>% as.data.frame () 
        out$arche.name <- unlist(out$arche.name)
        out <- out %>% utilsFanc::change.name.fanc("seqnames", "chr") %>% 
          utilsFanc::factor2character.fanc()
        return(out)
      }) %>% Reduce(rbind, .)
    aln.df <- aln.df %>% arrange(chr, start)
  }
  
  aln.df <- aln.df %>% filter(! arche.name %in% arche.exclude)
  utilsFanc::write.zip.fanc(df = aln.df, out.file = paste0(out.dir, "/", root.name, "_motif.tsv"), 
                            bed.shift = F, zip = F, col.names = T)
  
  if (to.bed == T) {
    # probably only works for consensus based scanning...
    aln.df.2 <- aln.df %>% mutate(chr="aln") %>% unique()
    trash <- aba.consensus.2.bed(abao = abao, aln.df = aln.df.2, out.dir = out.dir, root.name = root.name,
                                 split.by = split.by, split.by.int.name = split.by.int.name, ignore.strand = T)
  }
  return(aln.df)
}

aba.homer.scan <- function(abao, motif.files, bKeep.motif.seq = F,
                           out.dir,
                           findMotifs.pl = FINDMOTIFS.PL) {
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  if (is.null(names(motif.files))) {
    names(motif.files) <- basename(motif.files) %>% sub(".motif$", "", .)
  }
  motif.tmp.dir <- paste0(out.dir, "/motif_file_tmp/")
  dir.create(motif.tmp.dir, showWarnings = F)
  motif.mods <- lapply(names(motif.files), function(motif.name) {
    motif.file <- motif.files[motif.name]
    motif.mod <- paste0(motif.tmp.dir, "/", motif.name, ".motif")
    system(paste0("head -1 ", motif.file, " > ", motif.mod))
    header <- read.table(motif.mod, header = F, sep = "\t", quote = "")
    if (!bKeep.motif.seq) {
      header$V2 <- motif.name
    } else {
      header$V2 <- paste0(motif.name, ",", sub(",.+$", "", header$V2))
    }
    write.table(header, motif.mod, sep = "\t", quote = F, row.names = F, col.names = F)
    system(paste0("sed '1d' ", motif.file, " >> ", motif.mod))
    return(motif.mod)
  }) %>% unlist()
  
  motif.all <- paste0(motif.tmp.dir, "/all_motifs.motif")
  system(paste0("cat ", paste0(motif.mods, collapse = " "), " > ", motif.all))
  
  in.fa <- paste0(out.dir, "/", abao@meta.data$fa.root.name, "in.fa")
  seqinr::write.fasta(sequences = abao@ori.fa, names = names(abao@ori.fa),
                      file.out = in.fa)
  homer.out <- paste0(out.dir, "/", abao@meta.data$fa.root.name, "homer.tsv")
  cmd <- paste0(findMotifs.pl, " ", in.fa, " fasta ", out.dir, " -find ", motif.all, 
                " > ", homer.out)
  print(cmd)
  system(cmd)
  homer.df <- read.table(homer.out, sep = "\t", header = T, quote = "")
  colnames(homer.df) <- c("regionID", "offset", "sequence", "motif", "motif.strand", "score")
  homer.df <- inner_join(utilsFanc::gr2df(abao@ori.gr, keep.width = F), homer.df, by = "regionID") %>% 
    dplyr::mutate(center = ceiling((start + end) / 2)) %>% 
    dplyr::rename(region.strand = strand)
  
  homer.pos <- homer.df %>% dplyr::filter(motif.strand == "+") %>% 
    dplyr::mutate(start = center + offset, end = center + offset + nchar(sequence) - 1)
  homer.neg <- homer.df %>% dplyr::filter(motif.strand == "-") %>% 
    dplyr::mutate(end = center + offset, start = end - nchar(sequence) + 1)
  homer.df <- rbind(homer.pos, homer.neg)
  
  homer.df <- utilsFanc::df.rearrange.cols(df = homer.df, cols = "motif.strand", before = "region.strand")
  
  bed.df <- homer.df %>% dplyr::mutate(forth = paste0(sub(",.+$", "", motif), ";", sequence, ";", score),
                                       strand = motif.strand) %>% 
    dplyr::select(chr, start, end, forth, genome, strand)
  bed.df %>% split(., f = .$genome) %>% 
    lapply(function(df) {
      utilsFanc::write.zip.fanc(df, paste0(out.dir, "/", abao@meta.data$fa.root.name, "homer_",df$genome[1],".bed"), 
                                bed.shift = T)
    })
  
  parse.out <- paste0(out.dir, "/", abao@meta.data$fa.root.name, "homer_parse.tsv")
  write.table(homer.df, parse.out, row.names = F, col.names = T, sep = "\t", quote = F)
  return(homer.df)
}

aba.fimo.scan <- function(abao, motif.files, out.dir, thresh, fimo = FIMO) {
  # re-write aba.label.fimo to make it more like aba.homer.scan
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (length(motif.files) > 1) {
    motif.all <- paste0(out.dir, "/motif_all.meme")
    system(paste0("cat ", paste0(motif.files, collapse = " "), " > ", motif.all))
  } else {
    motif.all <- motif.files
  }
  in.fa <- paste0(out.dir, "/", abao@meta.data$fa.root.name, "in.fa")
  seqinr::write.fasta(sequences = abao@ori.fa, names = names(abao@ori.fa),
                      file.out = in.fa)
  fimo.out.dir <- paste0(out.dir, "/fimo")
  system("rm -rf " %>% paste0(fimo.out.dir))
  cmd <- paste0(fimo, " -o ",  fimo.out.dir, " --thresh ", thresh, " ",
                motif.all, " ", in.fa)
  print(cmd)
  system(cmd)
  fimo.out <- paste0(fimo.out.dir, "/fimo.tsv")
  fimo.df <- read.table(fimo.out, sep = "\t", header = T, quote = "")
  colnames(fimo.df) <- c("motif_id", "motif_alt_id", "regionID",
                         "start.on.region", "end.on.region", "motif.strand",
                         "score", "p.value", "q.value", "sequence")
  fimo.df <- inner_join(utilsFanc::gr2df(abao@ori.gr, keep.width = F), fimo.df, by = "regionID") %>% 
    dplyr::rename(region.strand = strand) %>% 
    dplyr::mutate(end = start + end.on.region - 1) %>% 
    dplyr::mutate(start = start + start.on.region - 1)
  fimo.df$log10p <- -1 * log10(fimo.df$p.value)
  fimo.df <- utilsFanc::df.rearrange.cols(df = fimo.df, cols = "motif.strand", before = "region.strand")
  
  bed.df <- fimo.df %>% dplyr::mutate(forth = paste0(motif_id, ";", motif_alt_id, ";", sequence, ";", score),
                                       strand = motif.strand) %>% 
    dplyr::select(chr, start, end, forth, genome, strand)
  bed.df %>% split(., f = .$genome) %>% 
    lapply(function(df) {
      utilsFanc::write.zip.fanc(df, paste0(out.dir, "/beds/", abao@meta.data$fa.root.name, "fimo_",df$genome[1],".bed"), 
                                bed.shift = T)
    })
  
  parse.out <- paste0(out.dir, "/", abao@meta.data$fa.root.name, "fimo_parse.tsv")
  write.table(fimo.df, parse.out, row.names = F, col.names = T, sep = "\t", quote = F)
  return(fimo.df)
}

aba.genomeAlign <- function(abao, ref.genome, query.genomes = NULL, ref.region = NULL,
                            bSLangan = F, SLangan.dir = "SLangan") {
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  if (is.null(query.genomes)) {
    query.genomes <- abao@ori.df$genome %>% .[. != ref.genome]
  }
  dups <- query.genomes[duplicated(query.genomes)]
  if (length(dups) > 0) {
    stop(paste0("some of the query genomes have more than one gene in the alignment: \n",
                paste0(dups, collapse = "\n")))
  }
  ori.df.ref <- abao@ori.df %>% filter(genome == ref.genome)
  if (nrow(ori.df.ref) < 1) {
    stop("nrow(ori.df.ref) < 1")
  }
  
  ori.df.query <- abao@ori.df %>% filter(genome %in% query.genomes)
  if (nrow(ori.df.query) < 1) {
    stop("nrow(ori.df.query) < 1")
  }
  
  gAln.dfs <- ori.df.query %>% split(., f = factor(.$genome, levels = .$genome)) %>% 
    lapply(function(q.region) {
      out.dir <- paste0(abao@work.dir, "/genomeAlign/")
      dir.create(out.dir, showWarnings = F, recursive = T)
      root.name <- paste0(abao@meta.data$fa.root.name, ref.genome, ".", q.region$genome)
      gAln.df <- ori.df.ref %>% split(., f = factor(.$regionID, levels = .$regionID)) %>% 
        lapply(function(r.region) {
          if (bSLangan) {
            aln.df <- aba.read.slangan.core(abao = abao, seq1 = r.region$regionID, seq2 = q.region$regionID,
                                            SLangan.dir = SLangan.dir)
          } else{
            aln.df <- abao@aln.fa[c(r.region$regionID, q.region$regionID)] %>% as.data.frame() 
          }
          aln.df <- aln.df %>% 
            filter(!!as.name(r.region$regionID) != "-" | !!as.name(q.region$regionID) != "-")

          bGap.mat <- aln.df != "-"
          # check if pairwise alignment is valid:
          lapply(list(r.region$regionID, q.region$regionID), function(rID) {
            len <- bGap.mat[, rID] %>% sum()
            len.ori <- abao@ori.df %>% dplyr::filter(regionID == rID) %>% 
              dplyr::mutate(width = end - start + 1) %>% 
              dplyr::pull(width)
            if (len != len.ori) {
              stop(paste0("invalide pairwise alignment for ", rID,
                          ": alignment length is ", len, 
                          " but length in abao@ori.df is ", len.ori))
            }
            return()
          })
          
          unGap.ids <- which(rowSums(bGap.mat) == ncol(bGap.mat))
          if (q.region$genome == "hg38" && r.region$regionID == "B6.Ly49c") {
            # browser()
          }
          gap.size.start <- unGap.ids[1] - 1 # size of gap at the start of pairwise alignment
          gap.size.end <- nrow(bGap.mat) - max(unGap.ids) # size of gap at the end of pairwise alignment
          
          bOut <- rowSums(bGap.mat) > 0 # should be all true
          if (gap.size.start > 0)
            bOut[1:gap.size.start] <- F # remove unaligned head and tail
          if (gap.size.end > 0)
            bOut[(length(bOut)-gap.size.end + 1):length(bOut)] <- F
          out.seq.df <- aln.df[bOut, ]
          # browser()
          
          grs <- lapply(list(r.region, q.region), makeGRangesFromDataFrame, keep.extra.columns = T)
          gr <- lapply(grs, function(gr) {
            if (gap.size.start < 1) {
              n.remove <- 0
            } else {
              n.remove <- sum(aln.df[, gr$regionID][1:gap.size.start] != "-")
            }
            gr <- resize(gr, width = width(gr) - n.remove, fix = "end", ignore.strand = F)
            
            if (gap.size.end < 1) {
              n.remove <- 0
            } else {
              n.remove <- sum(aln.df[, gr$regionID][(nrow(aln.df) - gap.size.end + 1):nrow(aln.df)] != "-")
            }
            gr <- resize(gr, width = width(gr) - n.remove, fix = "start", ignore.strand = F)
            return(gr)
          }) %>% Reduce(c, .)
          df.shifted <- gr %>% utilsFanc::gr2df(keep.width = F)
          
          cols <- c("chr", "start", "end", "regionID")
          gAln.df <- cbind(df.shifted[1, cols],
                           df.shifted[2, cols])
          colnames(gAln.df) <- c(cols, paste0("q.", cols))
          gAln.df$strand <- ifelse(df.shifted$strand[1] == df.shifted$strand[2], "+", "-")
          gAln.df$targetseq <- out.seq.df[, r.region$regionID] %>% paste0(collapse = "")
          gAln.df$queryseq <- out.seq.df[, q.region$regionID] %>% paste0(collapse = "")
          return(gAln.df)
        }) %>% Reduce(rbind, .)
      utilsFanc::write.zip.fanc(df = gAln.df, out.file = paste0(out.dir, "/", root.name, ".bed"), bed.shift = T)
      out.track <- paste0(out.dir, "/", root.name, ".align")
      gAln.df2track(gAln.df = gAln.df, out.track = out.track)
      out.track.gz <- paste0(out.track, ".gz")
      if (!is.null(ref.region)) {
        gAln.combine(gAln.1 = out.track.gz,
                     gAln.2 = paste0("~/genomes/", ref.genome, "/genomeAlign/", ref.genome,
                                     ".", q.region$genome, "/", ref.genome, ".", q.region$genome, 
                                     ".align.gz"),
                     ref.region = ref.region,
                     gAln.out = paste0(out.dir, "/", abao@meta.data$fa.root.name, "addTo_",
                                       ref.genome, ".", q.region$genome, ".align"))
      }
      
      return(gAln.df)
    }) 
  invisible(gAln.dfs)
}

gAln.df2track <- function(gAln.df, out.track = NULL) {
  out.df <- gAln.df[, 1:3]
  out.df$forth <- paste0("id:", 1:nrow(gAln.df),',genomealign:{chr:"', gAln.df$q.chr,
                         '",start:', gAln.df$q.start, ',stop:', gAln.df$q.end,
                         ',strand:"', gAln.df$strand, '",targetseq:"', gAln.df$targetseq,
                         '",queryseq:"', gAln.df$queryseq, '"}')
  if (!is.null(out.track)) {
    utilsFanc::write.zip.fanc(out.df, out.track, bed.shift = T)
  }
  return(out.df)
}

gAln.combine <- function(gAln.1, gAln.2, prefer.gAln.1 = T,
                         ref.region = NULL, gAln.out,
                         tabix = TABIX, keep.tmp = F) {
  if (!is.null(ref.region)) {
    tmp.1 <- paste0(gAln.out, "_tmp1")
    tmp.2 <- paste0(gAln.out, "_tmp2")
    system(paste0("tabix ", gAln.1, " ", ref.region, " > ", tmp.1))
    system(paste0("tabix ", gAln.2, " ", ref.region,  " > ", tmp.2))
    gAln.1 <- tmp.1
    gAln.2 <- tmp.2
  }
  gr1 <- rtracklayer::import.bed(gAln.1)
  gr2 <- rtracklayer::import.bed(gAln.2)
  o <- findOverlaps(query = gr1, subject = gr2, ignore.strand = T)
  if (prefer.gAln.1 && length(subjectHits(o)) > 0) {
    gr2 <- gr2[-subjectHits(o)]
  } else if (length(queryHits(o)) > 0) {
    gr1 <- gr1[-queryHits(o)]
  }
  gr <- c(gr1, gr2) %>% sort(ignore.strand = T)
  gr$name <- sub("^id[^,]+", "", gr$name)
  gr$name <- paste0("id:", 1:length(gr), gr$name)
  utilsFanc::write.zip.fanc(gr, gAln.out, bed.shift = T)
  if (!keep.tmp) {
    system(paste0("rm -rf ", gAln.out, "_tmp*"))
  }
  invisible(gr)
}
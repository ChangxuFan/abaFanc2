mini.genome.gen <- function(df, parent.genome, parent.genome.dir = "~/genomes/",
                            master.dir, threads = 4, 
                            bwa = BWA,
                            bowtie2.build = BOWTIE2.BUILD,
                            samtools = SAMTOOLS) {
  # df should be a bed format: chr, start, end, genome.name, parent.genome. parent.geome overrides global parameter
  # no headers
  
  # note: the header of the fasta file should be the chromosome name of the parent genome, so that it
  # can be shifted back to parent genome easily
  
  # parent.genome can be a fa file. In this case put parent.genome.dir to null
  
  if (is.character(df))
    df <- read.table(df, as.is = T, header = F)
  colnames(df) <- c("chr", "start", "end", "genome.name", "parent.genome")[1:ncol(df)]
  if (is.null(df$parent.genome)) {
    df$parent.genome <- parent.genome
  }
  if (is.null(df$genome.name)) {
    df$genome.name <- paste0(df$chr, "_", df$start, "_", df$end)
  }
  mini.genomes <- df %>% split(f = 1:nrow(df)) %>% 
    mclapply(function(x) {
      if (is.null(parent.genome.dir)) {
        parent.genome.fa <- x$parent.genome
      } else {
        parent.genome.fa <- paste0(parent.genome.dir, "/", parent.genome, "/", parent.genome, ".fa")
      }
      
      fa <- get.fasta.core(genome.fa = parent.genome.fa, chr = x$chr, 
                                    start = x$start, end = x$end)
      out.fa <- paste0(master.dir, "/", x$genome.name, "/", x$genome.name, ".fa")
      dir.create(dirname(out.fa), showWarnings = F, recursive = T)
      seqinr::write.fasta(sequences = fa, names = x$chr, file.out = out.fa)
      system(paste0(samtools, " faidx ", out.fa))
      system(paste0(bwa, " index ", out.fa))
      index <- liteRnaSeqFanc::bowtie2.index(fa = out.fa, threads = 1, bowtie2.build = bowtie2.build,
                                             log.file = paste0(dirname(out.fa), "/bowtie2-build.log"))
      trash <- utilsFanc::write.zip.fanc(df = x, out.file = paste0(master.dir, "/", x$genome.name, "/parent.bed"))
      return(x$genome.name)
    }, mc.cores = threads, mc.cleanup = T)
  return(mini.genomes)
}

mini.genomes.verify <- function(mini.genomes, master.dir, parent.genome, length.threshold,
                                parent.genome.fa = NULL, genome.dir = "~/genomes/",
                                blast.params = "", run.blast = T,
                                verify.dir) {
  # cat all the mini genomes into 1 fasta file
  dir.create(verify.dir, recursive = T, showWarnings = F)
  fa.vec <- paste0(master.dir, "/", mini.genomes, "/", mini.genomes, ".fa")
  fa.cat <- paste0(verify.dir, "/cat.fa")
  cmd <- paste0("cat ", paste0(fa.vec, collapse =  " "), " > ", fa.cat)
  system(cmd)
  ## now read in and change name:
  fa.cat.in <- seqinr::read.fasta(fa.cat, forceDNAtolower = F)
  seqinr::write.fasta(sequences = fa.cat.in, names = mini.genomes, file.out = fa.cat)
  rm(fa.cat.in)
  bed.vec <- paste0(master.dir, "/", mini.genomes, "/parent.bed")
  bed.cat <- paste0(verify.dir, "/cat.bed")
  cmd <- paste0("cat ", paste0(bed.vec, collapse =  " "), " > ", bed.cat)
  system(cmd)
  system(paste0("~/scripts/bed_browser_v2.sh ", bed.cat))
  cat(utilsFanc::bash2ftp(bed.cat))
  ### now align with blast
  if (is.null(parent.genome.fa))
    parent.genome.fa <- paste0(genome.dir, "/", parent.genome, "/", parent.genome, ".fa")
  blast.table <- paste0(verify.dir, "/blast.txt")
  
  if (run.blast == T) {
    trash <- abaFanc::blastn.fanc(query.fa.vec = fa.cat, subject.fa.vec = parent.genome.fa,
                                  outtable = blast.table, other.params = blast.params)
  }

  abaFanc::blast.table.to.bed(blast.table = blast.table, expand = F, squish = T, igv.gen = T,
                              push.2.igv = T, 
                              filter.func = function(x) return(x %>% abaFanc::blast.filter.length(threshold = length.threshold) %>% 
                                                                 abaFanc::blast.filter.best.match()))
  return()
}
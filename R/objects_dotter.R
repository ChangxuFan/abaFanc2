dotter.gff.from.bed <- function(fa.vec, bed.vec, out.gff) {
  # fa.vec: your x and y sequences. note it must follow the naming convention of
  # >name@chr:start-end
  # bed.vec: a vector of bed files, they will be read in and merged
  
  gr <- lapply(bed.vec, utilsFanc::import.bed.fanc, return.gr = T) %>% 
    Reduce(c, .)
  mcols(gr) <- NULL
  dir.create(dirname(out.gff), showWarnings = F, recursive = T)
  write("##gff-version 3", out.gff, sep = "\t")
  gtf.lines <- lapply(fa.vec, function(fa) {
    seq <- seqinr::read.fasta(fa)[[1]]
    loc.gr <- attr(seq, "Annot") %>% sub(".+ ", "", .) %>% 
      utilsFanc::loci.2.gr()
    include <- subsetByOverlaps(gr, loc.gr)
    include <- utilsFanc::gr2df(include)[, c("chr", "start", "end")]
    include$chr <- attr(seq, "name")
    # luckily, you don't need to consider strand issues.
    include$start <- pmax(0, include$start - start(loc.gr) + 1)
    include$end <- pmin(include$end - start(loc.gr) + 1, width(loc.gr))
    gtf.lines <- lapply(c("mRNA", "exon", "CDS"), function(type) {
      if (type == "mRNA") {
        use <- include[1, ]
        use$start <- min(include$start)
        use$end <- max(include$end)
      } else {
        use <- include
      }
      phase <- ifelse(type == "CDS", "0", ".")
      meta <- ifelse(type == "mRNA", paste0("ID=RNA", use$chr[1], ";Name=NAME", use$chr[1]),
                     paste0("Parent=RNA", use$chr[1]))
      df <- data.frame(chr = use$chr, source = "miao", type = type, 
                       start = use$start, end = use$end, score = ".", 
                       strand = "-", phase = phase, meta = meta)
    }) %>% do.call(rbind, .)
    header <- paste0("##sequence-region ", gtf.lines$chr[1], " ", 
                       min(gtf.lines$start), " ", max(gtf.lines$end))
    write(header, out.gff, sep = "\n", append = T)
    write.table(gtf.lines, out.gff, append = T, sep = "\t", quote = F, 
                col.names = F, row.names = F)
    return(gtf.lines)
  }) %>% do.call(rbind, .)
  invisible(gtf.lines)
}
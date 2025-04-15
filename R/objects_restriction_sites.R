fa.find.restriction.sites <- function(fa, RE.df, out.dir = NULL, root.name = NULL,
                                      ignore.revcomp = F,
                                      one.file = F) {
  # RE.df <- data.frame(enzyme = HindIII, seq = AAGCTT)
  utilsFanc::check.intersect(c("enzyme", "seq"), "required columns",
                             colnames(RE.df), "colnames(RE.df)")
  RE.df$seq <- toupper(RE.df$seq)
  RE.df$revcomp <- Biostrings::DNAStringSet(RE.df$seq) %>% Biostrings::reverseComplement() %>% 
    as.character() %>% toupper()
  
  utilsFanc::t.stat("Reading in fasta")
  
  seq <- seqinr::read.fasta(fa, forceDNAtolower = TRUE, as.string = T)
  
  if (is.null(out.dir)) out.dir <- dirname(fa)
  if (is.null(root.name)) root.name <- tools::file_path_sans_ext(basename(fa))
  
  pos <- lapply(1:nrow(RE.df), function(i) {
    enzyme <- RE.df[i, "enzyme"]
    utilsFanc::t.stat(">>>>>>>> Processing Enzyme: ", enzyme)
    
    motif <- RE.df[i, "seq"]
    motif.revcomp <- RE.df[i, "revcomp"]
    
    pos <- lapply(seq, function(seq) {
      chr <- attr(seq, "name")
      utilsFanc::t.stat("Processing: ", chr)
      seq <- as.character(seq) %>% toupper()
      
      pos <- stringr::str_locate_all(seq, motif)[[1]] %>% as.data.frame()
      
      if (motif != motif.revcomp && !ignore.revcomp) {
        stop("Finding revcomp part of the function is untested")
        pos$strand <- "+"
        
        pos.2 <- stringr::str_locate_all(seq, motif.revcomp)[[1]] %>% as.data.frame()
        pos.2$strand <- "-"
        
        pos <- rbind(pos, pos.2)
      }
      
      chr <- data.frame(chr = rep(chr, nrow(pos)))
      pos <- cbind(chr, pos)
      return(pos)
    }) %>% do.call(rbind, .)
    pos <- pos %>% arrange(chr, start, end)
    if (!one.file) {
      out.file <- paste0(out.dir, "/", root.name, "_", enzyme, "_", motif, ".bed")
      utilsFanc::write.zip.fanc(pos, out.file, bed.shift = T, ez = T)
    } else{
      pos$name <- enzyme
    }
    return(pos)
  }) %>% do.call(rbind, .) %>% dplyr::arrange(chr, start, end)
  
  if (one.file) {
    out.file <- paste0(out.dir, "/", root.name, "_sites", ".bed")
    utilsFanc::write.zip.fanc(pos, out.file, bed.shift = T, ez = T)
  }
  return()
}


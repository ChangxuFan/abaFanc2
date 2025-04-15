uniprot.get.peptides <- function(genes, lengths = c(9, 10), 
                                 uniprot.fa,
                                 return.strings = F) {
  utilsFanc::check.dups(genes, "genes")
  
  fa <- seqinr::read.fasta(uniprot.fa, forceDNAtolower = F)
  header <- names(fa)
  annot <- seqinr::getAnnot(fa) %>% unlist()

  gene.names <- stringr::str_extract(annot, "GN=[^ ]+ ") %>%
    gsub("GN=| +", "", .)
  
  utilsFanc::check.intersect(genes, "genes", gene.names, "genes in uniprot.fa")
  
  bUse <- gene.names %in% genes
  peptides <- fa.get.peptides(fa[bUse], lengths = lengths, return.strings = return.strings)
  names(peptides) <- gene.names[bUse]
  
  if (sum(duplicated(names(peptides))) > 0) {
    peptides <- peptides %>% split(f = names(peptides)) %>% 
      lapply(function(x) {
        if (length(x) > 1)
          return(Reduce(c, x))
        else
          return(x)
      })
  }
  if (!identical(sort(genes), sort(names(peptides)))) {
    stop("!identical(sort(genes), sort(names(peptides)))")
  }
  peptides <- peptides[genes]
  return(peptides)
}

fa.get.peptides <- function(fa, lengths, return.strings = F) {
  # fa <- list(seq1 = c("A", "K", "Y", "R", "D"),
  #            seq2 = c("D", "E", "W", "Y"))
  if (file.exists(fa[[1]][[1]])) {
    fa <- seqinr::read.fasta(fa, forceDNAtolower = F)
  }
  res <- lapply(fa, function(fa) {
    lapply(lengths, function(l.peptide) {
      l.fa <- length(fa)
      if (l.fa < l.peptide) {
        return()
      }
      
      last.pos <- l.fa - l.peptide + 1
      res <- lapply(1:last.pos, function(i) fa[i:(i + l.peptide - 1)])
      return(res)
    }) %>% Reduce(c, .) %>% return()
  })
  if (return.strings) {
    res <- lapply(res, function(each.gene) {
      unlist(lapply(each.gene, paste0, collapse = ""))
    })
  }
  return(res)
}

uniprot.header.genename.map <- function(headers, uniprot.fa) {
  fa <- seqinr::read.fasta(uniprot.fa, forceDNAtolower = F)
  headers.in.fa <- names(fa)
  
  utilsFanc::check.intersect(headers, "headers", headers.in.fa, "headers.in.fa")
  
  annot <- seqinr::getAnnot(fa) %>% unlist()
  gene.names <- stringr::str_extract(annot, "GN=[^ ]+ ") %>%
    gsub("GN=| +", "", .)
  
  if (length(gene.names) != length(headers.in.fa)) {
    stop("length(gene.names) != length(headers.in.fa)")
  }
  
  names(gene.names) <- headers.in.fa
  out <- gene.names[headers]
  return(out)
}

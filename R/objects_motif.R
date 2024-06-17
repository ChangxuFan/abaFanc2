homer.motif.2.pwMatrix <- function(motif.files) {
  # this is modeled after 
  # .requirePackage("chromVARmotifs", installInfo = "devtools::install_github(\"GreenleafLab/chromVARmotifs\")")

  motifs <- lapply(motif.files, function(motif.file) {
    lines <- readLines(motif.file)
    split <- rep(0, length(lines))
    split[grepl("^>", lines)] <- 1
    split <- cumsum(split)
    motifs <- lines %>% split(f = split)
    motifs <- lapply(motifs, function(motif) {
      title <- motif[1] %>% strsplit("\t") %>% unlist()
      title[2] <- sub(",BestGuess.+$", "", title[2])
      motif.matrix <- motif[-1] %>% strsplit("\t") %>% 
        lapply(function(x) {
          x <- as.numeric(x)
          x <- x/sum(x)
          }) %>% 
        do.call(rbind, .) %>% t() 
      motif.matrix <- log(motif.matrix/0.25)
      # based on https://github.com/GreenleafLab/chromVARmotifs/README.md
      rownames(motif.matrix) <- c("A", "C", "G", "T")
      TFBSTools::PWMatrix(ID = title[2], name = sub("/.+$", "", title[2]),
                          profileMatrix = motif.matrix)
    })
    return(motifs)
  }) %>% Reduce(c, .)
  names <- lapply(motifs, TFBSTools::ID)
  utilsFanc::check.dups(names, "motif IDs")
  names(motifs) <- names
  do.call(TFBSTools::PWMatrixList, motifs)
}

gr.add.motif <- function(gr, motif.dir = "~/motifs/homer/sth/each_motif/", motifs, 
                         return.gr = T, gr.extend = 0, threads = 1) {
  if (is.null(names(motifs))) {
    stop("motifs must be named. set names(motifs)")
  }
  if (gr.extend > 0) {
    gr.ext <- gr + gr.extend
  } else{
    gr.ext <- gr
  }
  motifs.df <- utilsFanc::safelapply(names(motifs), function(motif.name) {
    motif <- motifs[motif.name]
    motif.file <- paste0(motif.dir, "/", motif, ".bed")
    if (!file.exists(motif.file)) stop(paste0("motif file not found: ", motif.file))
    motif.gr <- rtracklayer::import(motif.file)
    o <- findOverlaps(gr.ext, motif.gr)
    res <- rep(0, length(gr))
    res[(1:length(gr)) %in% queryHits(o)] <- 1
    return(res)
  }, threads = threads) %>% as.data.frame()
  colnames(motifs.df) <- names(motifs)
  if(return.gr) {
    stop("return.gr not developped yet")
  } else {
    motifs.df <- cbind(data.frame(locus = utilsFanc::gr.get.loci(gr)), motifs.df)
  }
  
}
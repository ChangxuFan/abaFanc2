aba.mb.prep <- function(abao, type = "aligned") {
  if (is.character(abao)) {
    if (dir.exists(abao)) {
      abao <- paste0(abao, "/abao.Rds")
    }
    abao <- readRDS(abao)
  }
  if (type == "noGap") {
    fa <- paste0(abao@work.dir, "/", abao@meta.data$fa.root.name, "aligned.fa")
    fa <- fa.aln.remove.gap(fa)
  } else if (type == "aligned") {
    fa <- paste0(abao@work.dir, "/", abao@meta.data$fa.root.name, "aligned.fa")
  } else if (type == "shrink") {
    stop("not developed yet")
  } else {
    stop("type not found")
  }
  mb.dir <- paste0(abao@work.dir, "/mrbayes/")
  dir.create(mb.dir, showWarnings = F, recursive = T)
  out.nex <- paste0(mb.dir,  "/", basename(fa)) %>% sub(".fa$", ".nex", .)
  fa.2.nex(fa, out.file = out.nex, use.partition = F, add.mb.GTR = F)
  return(out.nex)
}
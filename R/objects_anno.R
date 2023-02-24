add.gff3 <- function(df, outdir) {
  system(paste0("mkdir -p ", outdir))
  df$gff3 <- paste0(outdir, "/", df$acc, ".gff3")
  for (i in 1:nrow(df)) {
    cmd <- paste0("wget -c ", "\"",
                  "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=",
                  df[i, "acc"], "\"", " -O ", df[i, "gff3"])
    print(cmd)
    code <- system(cmd, intern = F)
    if (code != 0)
      stop("download command not successful")
  }
  trash <- lapply(df$gff3, function(x) {
    if(!file.exists(x)) {
      stop(paste0(x, " was not successfully created"))
    }
  })
  return(df)
}


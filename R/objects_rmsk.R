rmsk.simple.gen <- function(rmsk.in, col.num, rmsk.out = NULL, zip = F,
                            filter.for = c("SINE", "LINE", "LTR", "DNA")) {
  rmsk <- read.table(rmsk.in, header = F, quote = "", sep = "\t")
  rmsk.simple <- rmsk[rmsk[, col.num] %in% filter.for,]
  
  if (is.null(rmsk.out)) {
    rmsk.out <- utilsFanc::insert.name.before.ext(name = rmsk.in, insert = "simple", delim = ".")
  }
  trash <- utilsFanc::write.zip.fanc(df = rmsk.simple, out.file = rmsk.out, zip = zip)
  return(rmsk.out)
}
gbk.parse.anderson <- function(gbk, chr.name, out.bed = NULL) {
  if (is.null(out.bed)) {
    out.bed <- paste0(gbk, ".bed")
  }
  gbk <- readChar(gbk, file.size(gbk))
  gbk <- genbankr::parseGenBank(text = gbk, ret.seq = F)
  gbk$FEATURES[[length(gbk$FEATURES)]]$note <- gbk$FEATURES[[length(gbk$FEATURES)]]$note %>% 
    sub("ORIGIN.+$", "", .)
  gr <- gbk$FEATURES %>% lapply(makeGRangesFromDataFrame, keep.extra.columns = T) %>% 
    Reduce(c, .) %>% `names<-`(NULL) %>% as.data.frame() %>%
    dplyr::mutate(seqnames = chr.name) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  gr <- gr[-1]
  mcols(gr) <- mcols(gr)[, c("note", "type")]
  utilsFanc::write.zip.fanc(gr, out.bed, bed.shift = T)
  return(out.bed)
}


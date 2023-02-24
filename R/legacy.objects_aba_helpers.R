# aba.plot.aln <- function(abao, tracks = NULL, p.in = NULL, x.lim = NULL, gap.only = F,
#                          smart.cut.df = NULL, plot = T, 
#                          plot.out = NULL,height=5, width=20, 
#                          atcg.size = NULL, 
#                          sub.align.out = NULL, ...) {
#   if (is.null(atcg.size))
#     atcg.size <- 5
#   if (!is.null(tracks)) {
#     aln.df <- lapply(tracks, function(track) {
#       regionID <- abao$track.df %>% filter(track.name == track) %>% pull(regionID)
#       aln.df <- abao$map[[regionID]] %>% mutate(y = track)
#       return(aln.df)
#     }) %>% Reduce(rbind,.)
#     y.levels <- tracks
#   } else {
#     aln.df <- lapply(names(abao$map), function(regionID) {
#       aln.df <- abao$map[[regionID]] %>% mutate(y = regionID)
#       return(aln.df)
#     }) %>% Reduce(rbind, .)
#     y.levels <- names(abao$map)
#   }
#   
#   if (!is.null(x.lim))
#     aln.df <- aln.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
#   if (!is.null(smart.cut.df)) {
#     cut.res <- aba.smart.cut(abao = abao, track.bp.df = aln.df, smart.cut.df = smart.cut.df, ...)
#     aln.df <- cut.res$out.df
#     aln.df$aln <- factor(aln.df$aln, levels = cut.res$levels)
#   }
#   
#   if (!is.null(sub.align.out)) {
#     # browser()
#     aln.prep <- paste0(sub.align.out, ".pre")
#     seqs <- aln.df %>% split(., f = factor(.$y, levels = y.levels)) %>% 
#       lapply(function(x) return(x$seq))
#     seqinr::write.fasta(sequences = seqs, names = names(seqs), file.out = aln.prep)
#     trash <- mafft.fanc(in.fa = aln.prep, aligned.fa = sub.align.out)
#   }
#   
#   if (plot == T) {
#     if (is.null(p.in))
#       p <- ggplot()
#     else
#       p <- p.in
#     
#     if (gap.only == T) {
#       aln.df <- aln.df %>% filter(seq %in% c("-", "N", "n"))
#       # browser()
#       p <- p + geom_tile(data = aln.df, 
#                          mapping = aes(x=aln, y = factor(y, levels = y.levels),
#                                        fill = seq),
#                          inherit.aes = F)
#     } else {
#       p <- p + geom_text(data = aln.df, 
#                          mapping = aes(x=aln, y = factor(y, levels = y.levels),
#                                        label = seq, color = seq),
#                          inherit.aes = F, size = atcg.size)
#     }
#     if (!is.numeric(aln.df$aln) && is.null(p.in)) {
#       vlines <- cut.res$cut.df.buffer$start
#       p <- p + scale_x_discrete(drop = F, breaks = factor(cut.res$breaks %>% c(vlines), levels = cut.res$levels)) +
#         geom_vline(xintercept = factor(vlines, levels = cut.res$levels))
#       
#     }
#     # if (!is.null(atcg.size)) {
#     #   p <- p + theme(text = element_text(size = atcg.size))
#     # }
#     if (!is.null(plot.out)) {
#       ggsave(plot.out, p, width = width, height = height, units = "in", device = "png", dpi = 100)
#     }
#     
#     return(p)
#   }
#   
#   return()
# }
# 
# import.bed.fanc <- function(bed, col.names = c("chr", "start", "end", "forth", "fifth", "strand"),
#                             return.gr = F, no.shift = F) {
#   if (grepl(".gz$", bed)) {
#     bed <- sub(".gz", "",bed)
#   }
#   df <- read.table(bed, as.is = T, sep = "\t", quote = "", header = F)
#   colnames(df) <- col.names[1:ncol(df)]
#   if (no.shift != F) {
#     df$start <- df$start + 1
#   }
#   if (return.gr == T) {
#     df <- GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = T)
#   }
#   return(df)
# }
# 
# aba.smart.cut <- function(abao, track.bp.df, smart.cut.df, n.breaks = 10) {
#   if (is.character(smart.cut.df))
#     smart.cut.df <- read.table(smart.cut.df, as.is = T, sep = "\t", quote = "", header = T)
#   cut.df <- aba.map.2.consensus(abao = abao, df = smart.cut.df)
#   cut.df.buffer <- cut.df %>% mutate(start = start - buffer.left, end = end + buffer.right) %>% 
#     select(start, end)
#   out.df <- cut.df.buffer %>% split(f= 1:nrow(cut.df.buffer)) %>% 
#     lapply(function(x) {
#       out.sub <- track.bp.df %>% filter(aln > x$start, aln < x$end)
#       return(out.sub)
#     }) %>% Reduce(rbind, .)
#   # browser()
#   levels <- lapply(1:nrow(cut.df.buffer), function(i) {
#     return(cut.df.buffer[i, "start"]:cut.df.buffer[i, "end"])
#   }) %>% Reduce(c, .) %>% sort()
#   n.levels <- length(levels)
#   break.interval <- ceiling(n.levels/n.breaks)
#   break.pos <- which(!duplicated(floor(1:n.levels/break.interval)))
#   breaks <- levels[break.pos]
#   return(list(out.df = out.df, cut.df = cut.df, cut.df.buffer = cut.df.buffer, levels = levels, breaks = breaks))
# }
# 
# aba.map.2.consensus <- function(abao, df) {
#   # expected format: chr, start, end, followed by other columns
#   if (is.character(df))
#     df <- read.table(df, as.is = T, sep = "\t", quote = "", header = T)
#   df$regionID <- gsub("[^A-Za-z0-9]", ".", df$regionID)
#   df <- df %>% split(f = 1:nrow(df)) %>% 
#     lapply(function(x) {
#       start.regionID <- abao$ori.df %>% filter(regionID == x$regionID) %>% pull(start)
#       start.on.region <- x$start - start.regionID + 1
#       end.on.region <- x$end - start.regionID + 1
#       x$start <- aba.map.core(abao = abao, regionID = x$regionID, pos = start.on.region)
#       x$end <- aba.map.core(abao = abao, regionID = x$regionID, pos = end.on.region)
#       return(x)
#     })  %>% Reduce(rbind, .)
#   df$chr <- "aln"
#   return(df)
# }
# 
# aba.map.core <- function(abao, regionID, pos) {
#   map <- abao$map[[regionID]]
#   # browser()
#   pos.out <- data.frame(ori = pos) %>% left_join(map) %>% pull(aln)
#   return(pos.out)
# }
# 
# aba.get.regionID <- function(abao, track) {
#   regionID <- abao$track.df %>% filter(track.name == track) %>% pull(regionID)
#   return(regionID)
# }
# 
# # aba.sub.aln <- function(abao, regionID = NULL, smart.cut.df) {
# #   print("miao")
# #   print("miao")
# # }
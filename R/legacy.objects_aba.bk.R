# aba.create <- function(aba.df, genome = NULL, fa = NULL, df.is.bed, 
#                        mafft = "/opt/apps/mafft/7.427/mafft", mafft.params = " --auto --reorder",
#                        work.dir, threads = 1) {
#   # fields for aba.df: chr	start	end	regionID	genome/fa  strand.
#   # strand seems okay to be missing.
#   system(paste0("mkdir -p ", work.dir))
#   if (is.character(aba.df)) {
#     if (df.is.bed == T) {
#       df.bed <- read.table(aba.df, sep = "\t", as.is = T)
#       aba.df <- df.bed[, 1:4]
#       colnames(aba.df) <- c("chr", "start", "end", "regionID")
#       if (!is.null(df.bed$V6))
#         aba.dfdf$strand <- df.bed$V6
#     } else {
#       aba.df <- read.table(aba.df, sep = "\t", as.is = T, quote = "", header = T)
#     }
#   }
#   
#   aba.df$regionID <- gsub("[^A-Za-z0-9]", ".", aba.df$regionID)
#   
#   if (!is.null(genome))
#     aba.df$genome <- genome
#   if (!is.null(fa))
#     aba.df$fa <- fa
#   ori.fa <- get.fasta.bed.2(df = aba.df, genome = genome, fa = fa, 
#                             df.is.bed = df.is.bed, threads = threads)
#   aln.fa <- mafft.fanc(in.fa = ori.fa, mafft = mafft, mafft.params = mafft.params,
#                        aligned.fa = paste0(work.dir, "/aligned.fa"))
#   map <- map.space(aligned.fa = aln.fa, threads = threads)
#   aba <- list(work.dir = work.dir, ori.df = aba.df, 
#               ori.gr = GenomicRanges::makeGRangesFromDataFrame(aba.df, keep.extra.columns = T),
#               ori.fa = ori.fa, aln.fa = aln.fa, map = map)
#   return(aba)
# }
# 
# map.space <- function(aligned.fa, threads =1) {
#   # if (is.character(in.fa))
#   #   in.fa <- read.fasta.fanc(in.fa)
#   if (is.character(aligned.fa))
#     aligned.fa <- read.fasta.fanc(aligned.fa)
#   
#   maps <- mclapply(aligned.fa, function(fa) {
#     # browser()
#     df <- data.frame(aln = 1:length(fa),
#                      seq = fa)
#     df2 <- df[df$seq != "-",]
#     df2$ori <- 1:nrow(df2)
#     df <- left_join(df, df2)
#     return(df)
#   }, mc.cores = threads, mc.cleanup = T)
#   return(maps)
# }
# 
# bp.df.gen <- function(data.gr, region.df=NULL) {
#   if (is.null(region.df))
#     region.df <- data.gr
#   if (!is.data.frame(region.df)) {
#     region.df <- region.df %>% `names<-`(NULL) %>% as.data.frame() 
#     region.df$chr <- region.df$seqnames
#   }
#     
#   start <- region.df$start %>% min()
#   end <- region.df$end %>% max()
#   
#   if (length(seqnames(data.gr) %>% unique()) > 1)
#     stop("only one chromosome is allowed")
#   if (length(region.df$chr %>% unique()) > 1)
#     stop("only one chromosome is allowed")
#   if (length(seqnames(data.gr)) > 0) {
#     if (as.character(region.df$chr[1]) != as.character(seqnames(data.gr)[1]))
#       stop("chromosome must match")
#   }
#   chr <- as.character(region.df$chr[1])
#   bp.gr <- data.frame(chr = chr, start = start:end, end = start:end) %>% 
#     GenomicRanges::makeGRangesFromDataFrame()
#   j.gr <- bp.gr %>% plyranges::join_overlap_left(data.gr)
#   return(j.gr)
# }
# 
# aba.add.track <- function(abao, track.df, prepend.regionID = F, threads = 1) {
#   # required fields: regionID	track.name	track.type	track.file
#   if (is.character(track.df))
#     track.df <- read.table(track.df, header = T, as.is = T, sep = "\t", quote="")
#   
#   track.df$regionID <- gsub("[^A-Za-z0-9]", ".", track.df$regionID)
#   track.df$track.name <- gsub("[^A-Za-z0-9]", ".", track.df$track.name)
#   if (prepend.regionID == T) {
#     track.df$track.name <- paste0(track.df$regionID, "..", track.df$track.name)
#   }
#   abao$track.df <- track.df
#   if (sum(duplicated(track.df$name)) > 0) {
#     stop("track name must be unique")
#   }
#   track.bp <- track.df %>% split(., f = factor(.$track.name, levels = .$track.name)) %>% 
#     mclapply(function(track) {
#       print(track$track.name)
#       region <- abao$ori.gr %>% plyranges::filter(regionID == track$regionID)
#       # browser()
#       if (tolower(track$track.type) %in% c("bdg", "bedgraph")) {
#         gr <- rtracklayer::import.bedGraph(con = track$track.file, which = region)
#         gr.bp <- bp.df.gen(data.gr = gr, region.df = region)
#         if (is.null(gr.bp$score))
#           gr.bp$score <- NA
#         gr.bp$score[is.na(gr.bp$score)] <- 0
#       } else if (tolower(track$track.type) %in% c("bed", "rmsk")) {
#         gr <- import.bed.fanc(bed = track$track.file, return.gr = T)
#         gr.bp <- bp.df.gen(data.gr = gr, region.df = region)
#         if (is.null(gr.bp$forth))
#           gr.bp$forth <- NA
#         if (is.null(gr.bp$fifth))
#           gr.bp$fifth <- NA
#       } else {
#         stop("only bedgraph, bed, and rmsk are supported right now")
#       }
#       return(gr.bp)
#     }, mc.cleanup = T, mc.cores = threads)
#   
#   abao$track.bp <- track.bp
#   track.bp.df <- mclapply(seq_along(track.bp), function(i) {
#     # track <- track.bp$Klra8..CAGE
#     track.name.i <- names(track.bp)[i]
#     track <- track.bp[[i]] %>% `names<-`(NULL) %>% mcols() %>% as.data.frame() 
#     track <- track %>% mutate(ori = 1:nrow(track), track.name = track.name.i) # %>% select(ori, score)
#     
#     regionID <- abao$track.df %>% filter(track.name == track.name.i) %>% pull(regionID)
#     track.type <- abao$track.df %>% filter(track.name == track.name.i) %>% pull(track.type)
#     map <- abao$map[[regionID]]
#     track.mapped <- left_join(track, map) # %>% pull(score)
#     if (tolower(track.type) %in% c("bdg", "bedgraph")) {
#       track.mapped <- track.mapped %>% filter(!is.na(score))
#     } else if (tolower(track.type) %in% c("bed", "rmsk")) {
#       track.mapped <- track.mapped %>% filter(!is.na(forth))
#     }
#     return(track.mapped)
#   }, mc.cores = threads, mc.cleanup = T) %>% Reduce(rbind, .)
#   # names(track.bp.df) <- names(track.bp)
#   # track.bp.df <- as.data.frame(track.bp.df)
#   #browser()
#   abao$track.bp.df <- track.bp.df
#   return(abao)
# }
# 
# get.dendro.order <- function(df, axis) {
#   if (axis == 2)
#     df <- t(df)
#   dend <- hclust(dist(df, diag = T))
#   order <- dend$labels[dend$order]
#   return(order)
# }
# 
# aba.plot.hm <- function(abao, track.type, use.mafft.order = F,
#                         abs = F, scale.row = F, cluster.tracks = F,
#                         remove.zero.tracks = F, add.nucleotide = F, gap.only = F, x.lim = NULL,
#                         n.per.row = NULL, plot.out = NULL, height=5, width=20,
#                         smart.cut.df = NULL) {
#   # arguments specific for bdg tracks: abs, scale.row, cluster.tracks, remove.zero.tracks
#   track.bp.df <- abao$track.bp.df
#   if (!is.null(x.lim)) {
#     track.bp.df <- track.bp.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
#   }
#   vline.df <- NULL
#   if (!is.null(smart.cut.df)) {
#     track.bp.df <- aba.smart.cut(track.bp.df = track.bp.df, smart.cut.df = smart.cut.df)
#   }
#   tracks.all <- track.bp.df %>% pull(track.name) %>% unique()
#   tracks <- tracks.all
#   # track.bp.df <- track.bp.df %>% filter(track.name %in% tracks)
#   if (use.mafft.order == T) {
#     tracks <- data.frame(regionID = names(abao$map)) %>% left_join(abao$track.df) %>% pull(track.name) %>% 
#       .[!is.na(.)]
#   }
#   if (track.type == "bdg") {
#     if (remove.zero.tracks == T) {
#       track.sum <- track.bp.df %>% group_by(track.name) %>% filter(!is.na(score)) %>% summarise(sum = sum(score))
#       tracks.non.zero <- track.sum %>% filter(sum != 0) %>% pull(track.name) %>% unique()
#       tracks <- tracks[tracks %in% tracks.non.zero]
#     } 
#     track.bp.df <- track.bp.df %>% filter(track.name %in% tracks)
#     if (abs == T) {
#       track.bp.df$score <- abs(track.bp.df$score)
#     }
#     if (scale.row == T) {
#       track.bp.df <- track.bp.df %>% group_by(track.name) %>% mutate(score = scale(score, center = F)) %>% 
#         ungroup() %>% as.data.frame()
#     }
#     
#     if (cluster.tracks == T) {
#       mat <- reshape2::acast(track.bp.df, formula = aln ~ track.name, value.var = "score")
#       tracks <- get.dendro.order(as.data.frame(mat), axis = 2)
#     }
#     
#     p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = score)) +
#       geom_tile() +
#       scale_fill_gradient(low="white", high="red") 
#   }
#   
#   if (track.type %in% c("bed", "rmsk")) {
#     if (track.type == "bed") {
#       p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = forth)) +
#         geom_tile() 
#     }
# 
#     if (track.type == "rmsk") {
#       track.bp.df$fifth <- sub("/.+$", "", track.bp.df$fifth)
#       p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = fifth)) +
#         geom_tile() +
#         scale_fill_manual(values = c("darkorange", "grey20", "seagreen4", "grey30", "grey40", "red", "grey50", "blue"),
#                           breaks = c("LINE", "Low_complexity", "LTR", "Satellite", "Simple_repeat", "SINE", "Unknown", "DNA"))
#     }
#   }
#   
#   if (add.nucleotide == T) {
#     p <- aba.plot.aln(abao, tracks = tracks, p = p, x.lim = x.lim, gap.only = gap.only)
#   }
#   if (!is.null(plot.out)) {
#     ggsave(plot.out, p, width = width, height = height, units = "in", device = "png", dpi = 100)
#   }
#   return(p)
# 
# }
# 
# aba.plot.aln <- function(abao, tracks = NULL, p = NULL, x.lim = NULL, gap.only = F,
#                          plot.out = NULL,height=5, width=20) {
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
# 
#   if (is.null(p))
#     p <- ggplot()
#   
#   if (gap.only == T) {
#     aln.df <- aln.df %>% filter(seq %in% c("-", "N", "n"))
#     # browser()
#     p <- p + geom_tile(data = aln.df, 
#                        mapping = aes(x=aln, y = factor(y, levels = y.levels),
#                                      fill = seq),
#                        inherit.aes = F)
#   } else {
#     p <- p + geom_text(data = aln.df, 
#                        mapping = aes(x=aln, y = factor(y, levels = y.levels),
#                                      label = seq, color = seq),
#                        inherit.aes = F)
#   }
#   
#   if (!is.null(plot.out)) {
#     ggsave(plot.out, p, width = width, height = height, units = "in", device = "png", dpi = 100)
#   }
# 
#   return(p)
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
# aba.smart.cut <- function(track.bp.df, smart.cut.df) {
#   if (is.character(smart.cut.df))
#     smart.cut.df <- read.table(smart.cut.df, as.is = T, sep = "\t", quote = "", header = T)
#   cut.df <- aba.map.2.consensus(abao = abao, smart.cut.df = smart.cut.df)
#   cut.df.buffer <- cut.df %>% mutate(start = start - buffer.left, end = end + buffer.right) %>% 
#     select(start, end)
#   out.df <- cut.df %>% split(f= 1:nrow(cut.df.buffer)) %>% 
#     lapply(function(x) {
#       out.sub <- track.bp.df %>% filter(aln > x$start, aln < x$end)
#       return(out.sub)
#     }) %>% Reduce(rbind, .)
#   return(list(out.df = out.df, cut.df = cut.df))
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

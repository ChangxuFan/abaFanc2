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
#         aba.df$strand <- df.bed$V6
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
# aba.add.track <- function(abao, track.df, prepend.regionID = F, threads = 1,
#                           track.type) {
#   # required fields: regionID	track.name	track.type	track.file
#   if (is.character(track.df))
#     track.df <- read.table(track.df, header = T, as.is = T, sep = "\t", quote="")
#
#   track.df$regionID <- gsub("[^A-Za-z0-9]", ".", track.df$regionID)
#   track.df$track.name <- gsub("[^A-Za-z0-9]", ".", track.df$track.name)
#   track.df <- track.df %>% filter(regionID %in% abao$ori.df$regionID)
#
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
#       if (track.type == "bdg") {
#         gr <- rtracklayer::import.bedGraph(con = track$track.file, which = region)
#         gr.bp <- bp.df.gen(data.gr = gr, region.df = region)
#         if (is.null(gr.bp$score))
#           gr.bp$score <- NA
#         gr.bp$score[is.na(gr.bp$score)] <- 0
#       } else if (track.type %in% c("bed", "rmsk")) {
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
#   abao[[paste0("track.bp.", track.type)]] <- track.bp
#   track.bp.df <- mclapply(seq_along(track.bp), function(i) {
#     # track <- track.bp$Klra8..CAGE
#     track.name.i <- names(track.bp)[i]
#     track <- track.bp[[i]] %>% `names<-`(NULL) %>% mcols() %>% as.data.frame()
#     track <- track %>% mutate(ori = 1:nrow(track), track.name = track.name.i) # %>% select(ori, score)
#
#     regionID <- abao$track.df %>% filter(track.name == track.name.i) %>% pull(regionID)
#     map <- abao$map[[regionID]]
#     track.mapped <- left_join(track, map) # %>% pull(score)
#     if (track.type == "bdg") {
#       track.mapped <- track.mapped %>% filter(!is.na(score))
#     } else if (track.type %in% c("bed", "rmsk")) {
#       track.mapped <- track.mapped %>% filter(!is.na(forth))
#     }
#     return(track.mapped)
#   }, mc.cores = threads, mc.cleanup = T) %>% Reduce(rbind, .)
#   # names(track.bp.df) <- names(track.bp)
#   # track.bp.df <- as.data.frame(track.bp.df)
#   #browser()
#   abao[[paste0("track.bp.df.",track.type)]]  <- track.bp.df
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
#                         broadcast = F, fill.gap = F,
#                         abs = F, scale.row = F, cluster.tracks = F,
#                         remove.zero.tracks = F, x.lim = NULL,
#                         plot.out = NULL, height=5, width=20, text.size = NULL,
#                         add.nucleotide = F, gap.only = F,
#                         atcg.size = NULL, sub.align.out = NULL,
#                         smart.cut.df = NULL) {
#   # arguments specific for bdg tracks: abs, scale.row, cluster.tracks, remove.zero.tracks
#   track.bp.df <- abao[[paste0("track.bp.df.", track.type)]]
#   if (!is.null(x.lim)) {
#     track.bp.df <- track.bp.df %>% filter(aln >= x.lim[1], aln <= x.lim[2])
#   }
#
#   if (!is.null(smart.cut.df)) {
#     cut.res <- aba.smart.cut(abao = abao, track.bp.df = track.bp.df, smart.cut.df = smart.cut.df)
#     track.bp.df <- cut.res[[1]]
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
#
#     if (cluster.tracks == T) {
#       non.all.zero <- track.bp.df %>% na.omit() %>% group_by(aln) %>%
#         summarise(nz = sum(score != 0) > 0) %>% filter(nz == T) %>% pull(aln)
#
#       track.scale.center <- track.bp.df %>% filter(aln %in% non.all.zero) %>%
#         group_by(track.name) %>% mutate(score = scale(score, center = T)) %>%
#         ungroup() %>% as.data.frame()
#       mat <- reshape2::acast(track.scale.center, formula = aln ~ track.name, value.var = "score")
#       tracks <- get.dendro.order(as.data.frame(mat), axis = 2)
#     }
#
#     if (scale.row == T) {
#       track.bp.df <- track.bp.df %>% group_by(track.name) %>% mutate(score = scale(score, center = F)) %>%
#         ungroup() %>% as.data.frame()
#     }
#
#     if (!is.null(smart.cut.df)) {
#       track.bp.df$aln <- factor(track.bp.df$aln, levels = cut.res$levels)
#     }
#     p <- ggplot(track.bp.df, aes(x = aln, y = factor(track.name, levels = tracks), fill = score)) +
#       geom_tile() +
#       scale_fill_gradient(low="white", high="red")
#     # browser()
#   }
#
#   if (track.type %in% c("bed", "rmsk")) {
#     if (fill.gap == T) {
#       #browser()
#       track.bp.df <- track.bp.df %>% split(., f= factor(.$track.name, levels = unique(.$track.name))) %>%
#         lapply(function(df) {
#           regionID <- aba.get.regionID(abao, df$track.name[1])
#           feature.bound <- df %>% group_by(forth) %>% summarise(start = min(aln), end = max(aln)) %>%
#             ungroup() %>% as.data.frame()
#           j <- left_join(abao$map[[regionID]], df)
#           filled <- feature.bound %>% split(., f=1:nrow(.)) %>%
#             lapply(function(feature) {
#               filled.sub <- j %>% filter(aln >= feature$start, aln <= feature$end)
#               filled.sub$forth <- feature$forth
#               filled.sub$track.name <- df$track.name[1]
#               return(filled.sub)
#             }) %>% Reduce(rbind, .)
#           return(filled)
#         }) %>% Reduce(rbind,.)
#     }
#
#     if (broadcast == T) {
#       track.bp.df$track.name <- "consensus"
#       track.bp.df <- track.bp.df %>% unique()
#       tracks <- "consensus"
#     }
#
#     if (!is.null(smart.cut.df)) {
#       track.bp.df$aln <- factor(track.bp.df$aln, levels = cut.res$levels)
#     }
#
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
#   if (!is.numeric(track.bp.df$aln)) {
#     vlines <- cut.res$cut.df.buffer$start
#     p <- p + scale_x_discrete(drop = F, breaks = factor(cut.res$breaks %>% c(vlines), levels = cut.res$levels)) +
#       geom_vline(xintercept = factor(vlines, levels = cut.res$levels))
#
#   }
#   if (!is.null(text.size)) {
#     p <- p + theme(text = element_text(size = text.size))
#   }
#
#   if (add.nucleotide == T) {
#     # browser()
#     p <- aba.plot.aln(abao, tracks = tracks, p.in = p, x.lim = x.lim, gap.only = gap.only,
#                       smart.cut.df = smart.cut.df, atcg.size = atcg.size, sub.align.out = sub.align.out)
#   }
#   if (!is.null(plot.out)) {
#     ggsave(plot.out, p, width = width, height = height, units = "in", dpi = 100)
#   }
#   return(p)
#
# }
#

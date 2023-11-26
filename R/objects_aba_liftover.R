aba.liftover <- function(abao, bed, is.bw = F, regionID.in, regionID.out = NULL, regions.exclude = NULL,
                         output.each = F, out.file = NULL) {
  # this is just a linear mapping, can't flip strands...
  # regionID.in must be scalar. regionID.out can be a vector
  abao <- aba.import(abao)
  in.chr <- abao@ori.df[abao@ori.df$regionID == regionID.in, "chr"]

  if (is.character(bed)) {
    if (is.bw) {
      bed <- rtracklayer::import(bed) 
      bed <- subsetByOverlaps(bed, abao@ori.gr[abao@ori.gr$regionID == regionID.in])
      bed <- bed %>% utilsFanc::gr2df()
    } else {
      bed <- import.bed.fanc(bed = bed, return.gr = F, no.shift = F)[, 1:6]
    }
  }

  if (is.bw) {
    bed$width <- NULL
    bed$strand <- NULL
    bed$forth <- utilsFanc::gr.get.loci(bed)
  }

  bed <- bed[bed$chr == in.chr,]
  if (nrow(bed) < 1) {
    stop("nrow(bed) < 1")
  }

  df.on.aln <- aba.map.bed.2.consensus(abao = abao, bed = bed, regionID = regionID.in)
  df.lifted <- aba.map.consensus.2.bed(abao = abao, bed = df.on.aln, regions.out = regionID.out, 
                                       regions.exclude = regions.exclude,
                                       output.each = output.each, out.file = out.file)
  return(df.lifted)
}

aba.map.bed.2.consensus <- function(abao, bed, regionID) {
  # this function doesn't deal with strand. It preserves strand info as is.
  if (is.character(bed)) {
    bed <- import.bed.fanc(bed = bed, return.gr = F, no.shift = F)[, 1:6]
  }
  bed.aln <- bed %>% split(., f= 1:nrow(.)) %>% 
    lapply(function(x) {
      ori.df <- abao@ori.df %>% .[.$regionID == regionID, ]
      
      if (x$chr != ori.df$chr)
        stop("chromosome doesn't match for " %>% paste0(x$forth))
      
      x$start <- max(ori.df$start, x$start)
      x$end <- min(x$end, ori.df$end)
      
      map <- abao@map[[regionID]]
      map <- map %>% dplyr::mutate(ori.genome = ori + ori.df$start - 1) %>% 
        .[!is.na(.$ori.genome),]
      for (y in c("start", "end")) {
        x[, y] <- map %>% .[.$ori.genome == x[, y], "aln"]
      }
      return(x)
    }) %>% Reduce(rbind, .)

  return(bed.aln)
}

aba.map.consensus.2.bed <- function(abao, bed, regions.out = NULL, regions.exclude = NULL,
                                    output.each = F, out.file = NULL) {
  # bed: bed format, but relative to the coordinates of aln
  if (is.null(regions.out)) {
    regions.out <- abao@ori.df$regionID %>% as.character()
  }
  if (!is.null(regions.exclude)) {
    regions.out <- regions.out %>% .[!grepl(paste0(regions.exclude, collapse = "|"),.)]
  }
  if (is.character(bed)) {
    bed <- import.bed.fanc(bed = bed, return.gr = F, no.shift = F)
  }
  if (! "forth" %in% colnames(bed)) {
    stop('! "forth" %in% colnames(bed)')
  }
  df <- bed %>% split(., f= 1:nrow(.)) %>% 
    lapply(function(x) {
      df <- lapply(regions.out, function(regionID) {
        ori.df <- abao@ori.df %>% .[.$regionID == regionID, ]
        if (x$chr != ori.df$chr) {
          warning("chromosome doesn't match for " %>% paste0(x$forth, " and region ", regionID))
          # return()
        }
        map <- abao@map[[regionID]]
        map <- map %>% dplyr::mutate(ori.genome = ori + ori.df$start - 1) %>% 
          .[!is.na(.$ori.genome),]
        
        inRegion <- map %>% .[.$aln >= x[, "start"] & .$aln <= x[, "end"], "ori.genome"]
        if (length(inRegion) == 0) {
          print(paste0("Not mapped: aln coordinates start: ", x$start, " ; end: ", x$end))
          return()
        } else {
          x$start <- min(inRegion)
          x$end <- max(inRegion)
        }
        
        x$chr <- ori.df$chr
        return(x)
      }) %>% Reduce(rbind, .)
      if (!is.null(out.file) && output.each == T) {
        id <- x$forth
        out.file <- out.file %>% sub(".bed$", paste0("_", id, ".bed"), .)
        utilsFanc::write.zip.fanc(df, out.file, bed.shift = T)
      }
      return(df)
    }) %>% Reduce(rbind, .)
  if (!is.null(out.file)) {
    utilsFanc::write.zip.fanc(df, out.file, bed.shift = T)
  }
  return(df)
}

aba.map.core <- function(abao, regionID, pos = NULL, aln = NULL) {
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  map <- abao@map[[regionID]]
  # browser()
  if (!is.null(pos))
    pos.out <- data.frame(ori = pos) %>% left_join(map) %>% pull(aln)
  else 
    pos.out <- data.frame(aln = aln) %>% left_join(map) %>% pull(ori)
  return(pos.out)
}



aba.shrink.2.full <- function(abao, pos, cons.name) {
  shrink.map <- abao@meta.data$consensus[[cons.name]]$shrink.map
  if (is.null(shrink.map))
    stop("is.null(shrink.map)")
  df <- data.frame(pos.shrink = pos)
  df <- left_join(df, shrink.map)
  pos.out <- df$pos.full
  if (length(pos) != length(pos.out))
    stop("length(pos) != length(pos.out)")
  return(pos.out)
}

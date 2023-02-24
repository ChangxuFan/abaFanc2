# the following function aims to get all Ly49 sequences from exon4 to pro1. 
# the basic strategy: blast exon4 and pro1. Then go through the blast result to get intervals
# The principle is: exon4 is the anchor. Each exon4 will correspond to 1 final interval. Pro1
#is the target, only those satisfying specified conditions will be used.
#for orphan exon4 with no paired pro1, just give it a random fixed-length interval.
range.gen.core <- function(beds.left, beds.right, direction = "+", dist,
                           search.start = 0, search.stop = NULL,
                           left.buffer=0, right.buffer=0,
                           out.bed = NULL, add.add.cols = NULL) {
  bed.left <- lapply(beds.left, function(bed.left) {
    if (is.character(bed.left))
      bed.left <- read.table(bed.left, as.is = T)
    return(bed.left)
  }) %>% Reduce(rbind, .)

  bed.right <- lapply(beds.right, function(bed.right) {
    if (is.character(bed.right))
      bed.right <- read.table(bed.right, as.is = T)
    return(bed.right)
  }) %>% Reduce(rbind, .)
  
  vec.left <- bed.left[, 2] %>% sort()
  vec.right <- bed.right[, 2] %>% sort()
  chr <- bed.left[1,1]
  
  if (direction == "+") {
    sign <- 1
    anchors <- vec.left
    targets <- vec.right
  } else {
    stop("this function has not been tested to work with direction = '-'")
    sign <- -1
    anchors <- vec.right
    targets <- vec.left
  }

  
  if (is.null(search.stop)) 
    search.stop <- anchors[length(anchors)]
  
  anchors <- anchors %>% .[. >= search.start & . <= search.stop]
  
  df <- lapply(anchors, function(anchor){
    range <- c(anchor, anchor + sign*dist) %>% sort()
    targets.in.range <- targets %>% .[. >= range[1] & . <= range[2]]
    if (length(targets.in.range) < 1) {
      target.found <- anchor + sign*dist
    } else {
      anchor.target.dists <- abs(targets.in.range - anchor)
      target.found <- targets.in.range[which.min(anchor.target.dists)]
    }
    interval <- c(anchor, target.found) %>% sort()
    interval[1] <- interval[1] - left.buffer
    interval[2] <- interval[2] + right.buffer
    df <- data.frame(chr = chr, start = interval[1], end = interval[2])
    return(df)
  }) %>% Reduce(rbind, .)
  if (!is.null(add.add.cols)) {
    df <- cbind(df, add.add.cols)
  }

  
  if (!is.null(out.bed))
    utilsFanc::write.zip.fanc(df = df, out.file = out.bed, zip = T)
  return(df)
}
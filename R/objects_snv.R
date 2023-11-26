tsv.2.snv <- function(tsv.file, snv, ref.name, shift) {
	# mostly rewriting ~/viralBrowser/release_github_4_22/tsv2snv.R
	tsv <- readr::read_tsv(tsv.file, col_types = "ccc", na = character())
	tsv$type <- apply(tsv, 1, function(x) {
	  types <- c("mismatch", "insertion", "deletion")
	  type <- types[which(x!="")]
	  if (length(type) > 1) {
	    if (!(all(type == c("mismatch", "insertion"))) )
	      stop (paste0(paste0(type, collapse = "-")," is a weird combo"))
	    else type <- "insertion"
	  }

	  if (length(type) ==0)
	    type <- ""
	  return(type)
	})

	tsv$chr <- rep(ref.name, nrow(tsv))
	tsv$left <- (1:nrow(tsv)) -1 + shift
	tsv$right <- tsv$left + 1

	tsv.out <- tsv %>% filter(type != "") %>% mutate(var = paste0(type ,": ",insertion, mismatch, deletion)) %>%
	  dplyr::select(chr, left, right, var)
	dir.create(dirname(snv), showWarnings = F, recursive = T)
	utilsFanc::write.zip.fanc(tsv.out, snv)
	return()
}
abao.2.snv <- function(abao, ref, query, keep.tempt = F, write.ori = F) {
	source("~/viralBrowser/v2/R/align.R")
	source("~/viralBrowser/v2/values_WangLabServer.R")
	source("~/viralBrowser/v2/R/pipe_ctrl.R")
	if (length(ref) != 1) {
		stop("length(ref) != 1")
	}
	
	utilsFanc::check.intersect(c(ref,query), "c(ref, query)", 
			names(abao@ori.fa), "names(abao@ori.fa)")
	
	out.dir <- paste0(abao@work.dir, "/snv/")
	
	tempt.dir <- paste0(out.dir, "/tmp/")
	dir.create(tempt.dir, showWarnings = F, recursive = T)
	ref.fa <- paste0(tempt.dir, "/ref.fa")
	seqinr::write.fasta(sequences = abao@ori.fa[ref], names = ref, file.out = ref.fa)

	lapply(query, function(query) {
		snv <- paste0(out.dir, "/ref_", ref, "..q_", query, ".bed")
		tempt.dir <- paste0(tempt.dir, "/", query, "/")
		dir.create(tempt.dir, showWarnings = F, recursive = T)
		
		query.fa <- paste0(tempt.dir, "/query.fa")
		seqinr::write.fasta(sequences = abao@ori.fa[query], names = query, file.out = query.fa)
		aln.fa <- paste0(tempt.dir, "/aln.fa")
		tsv <- paste0(tempt.dir, "/aln.tsv")

		align.stretcher(in.fa = query.fa, ref.fa = ref.fa, aln.fa = aln.fa)
		align.parse.markx3(aln.fa = aln.fa, out.tsv = tsv, write.ori = write.ori)
		tsv.2.snv(tsv.file = tsv, snv = snv, 
			ref.name = abao@ori.df[abao@ori.df$regionID == ref, "chr"],
		 	shift = abao@ori.df[abao@ori.df$regionID == ref, "start"] - 1)
		})
	
	if (!keep.tempt) {
		system(paste0("rm -rf ", tempt.dir))
	}
}

snv.2.abao <- function(snv, genome, buffer.1side = 3000, work.dir) {
	qname <- basename(snv) %>% sub(".gz", "", .) %>%
		tools::file_path_sans_ext()
	
	snv.df <- utilsFanc::import.bed.fanc(snv)
	start <- min(snv.df$start)
	end <- max(snv.df$end)
	if (length(unique(snv.df$chr)) != 1) {
		stop("only 1 chromosome is allowed in the snv file!")
	}
	chr <- snv.df$chr[1]
	if (end - start > 50000) {
		stop("region larger than 50000. Too large!")
	}

	Start <- start - buffer.1side
	End <- end + buffer.1side

	# r means the coordinates on reference, r means the coordinates on query
	region.r <- data.frame(chr = chr, start = Start, end = End, qname = qname)
	snv.q <- snv.df %>% dplyr::mutate(start = start - Start + 1, end = end - Start + 1)
	snv.region.q <- data.frame(chr = qname, start = min(snv.q$start), end = max(snv.q$end))

	fa <- get.fasta.gr(gr = makeGRangesFromDataFrame(region.r, keep.extra.columns = T),
		id.col = "qname", genome = genome, 
		drop.strand = T, wrap.fa = T, return.fa = T, as.string = F)
	# this is a seqinr object
	fa <- as.character(fa[[1]])
	
	map <- data.frame(ori = 1:length(fa), seq = fa)
	
	for (i in 1:nrow(snv.q)) {
		df <- snv.q[i, ]
		
		if (grepl("mismatch", df$forth)) {
			map[map$ori == df$start, "seq"] <- stringr::str_extract(df$forth, ".$")
		} else if (grepl("insertion", df$forth)) {
			ins <- df$forth %>% sub("insertion: ", "", .) %>% strsplit("") %>% unlist()
			
			map <- rbind(map[map$ori <= df$start,], data.frame(ori = df$start, seq = ins),
				map[map$ori > df$start,])

			rm(ins)
		} else if (grepl("deletion", df$forth)) {
			del <- df$forth %>% sub("deletion: ", "", .) %>% strsplit("") %>% unlist()
			map <- rbind(map[map$ori < df$start, ], 
						 map[map$ori >= df$start + length(del), ])
			rm(del)
		}
	}
	rm(i)
	rm(df)
	
	fa <- map$seq
	q.fa <- paste0(work.dir, "/", qname, ".fa")
	system(paste0("rm -rf ", work.dir))
	dir.create(work.dir, showWarnings = F, recursive = T)
	seqinr::write.fasta(sequences = list(fa), names = qname, file.out = q.fa)
	region.q <- data.frame(chr = qname, start = 1, end = length(fa), regionID = qname, genome = NA, fa = q.fa)
	region.r <- cbind(region.r[, c("chr", "start", "end")], 
		data.frame(regionID = "ref", genome = genome, fa = NA))
	
	np <- snv.region.q[, c("chr", "start", "end")] %>% 
		dplyr::mutate(V4 = ".", V5 = ".", V6 = ".", V7 = ".", V8 = ".", V9 = ".", 
			V10 = floor((end - start)/2))
	write.table(np, paste0(work.dir, "/", qname, ".narrowPeak"), sep = "\t", 
		col.names = F, row.names = F, quote = F)
	try({aba.create(aba.df = rbind(region.r, region.q), df.is.bed = F, work.dir = work.dir)})
	try(system(paste0("cut -f1-2 ", work.dir, "/", qname, ".fa.fai",
						" > ", work.dir, "/", qname, ".chrom.sizes")))
	return()
}



phyml.GTR.GIF <- function(in.phy, 
                      phyml = "/bar/cfan/software/PhyML-3.1/phyml", run = T) {
  if (grepl(".f(n)*a(sta)*$", in.phy)) {
    in.phy <- fasta2phylip(in.phy)
  }
  out.dir <- paste0(dirname(in.phy), "/phyml/")
  softlink <- paste0(out.dir, "/", basename(in.phy))
  dir.create(out.dir, showWarnings = F, recursive = T)
  system(paste0("ln -s ", normalizePath(in.phy), " ", softlink))
  
  cmd <- paste0(phyml, " -i ", softlink, " -d nt -o tlr -m GTR -c 4 -a e -v e -f e ",
                " -s SPR --rand_start --n_rand_starts 6 -b -4 ",
                " --quiet")
  cat(cmd)
  cat("\n")
  if (run)
    system(cmd)
  return(cmd)
}

phyml.titrate <- function(in.file, params.const, params.titrate, 
                          threads = NULL, base.params = "-d nt -o tlr --quiet -f e ",
                          phyml = PHYML) {
  out.dir <- paste0(dirname(in.file), "/phyml_",
                    paste0(params.const, collapse = "") %>% gsub(" +", "", .),
                    "/")
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (is.null(threads)) {
    threads <- min(threads, 16)
  }
  utilsFanc::safelapply(params.titrate, function(params) {
    softlink <- paste0(out.dir, "/", basename(in.file), "_",
                       paste0(params, collapse = "") %>% gsub(" +", "", .))
    system(paste0("ln -s ", normalizePath(in.file), " ", softlink))
    cmd <- paste0(phyml, " -i ", softlink, " ", base.params, " ",
                  paste0(params.const, collapse = " "), " ", 
                  paste0(params, collapse = " ")) %>% gsub(" +", " ",.)
    print(cmd)
    system(cmd)
  }, threads = threads)
  return()
}

phyml.build.models <- function() {
  models <- c("GTR", "HKY85", "K80", "TN93") %>% paste0("-m ",.)
  I <- c("", "-v e")
  G <- c("", "-a e -c 4")
  string <- outer(models, I, paste, sep = " ") %>% as.vector()
  string <- outer(string, G, paste, sep = " ") %>% as.vector()
  return(string)
}

phyml.build.trees <- function(n.bootstrap = NULL, branch.supports = c(1,2,4,5)) {
  s <- c("-s SPR --rand_start --n_rand_starts 5", "-s NNI")
  b <- c()
  if (!is.null(branch.supports)) {
    b <- c(b, paste0("-b ", -1 * branch.supports))
  }
  
  if (!is.null(n.bootstrap)) {
    b <- c(b, paste0("-b ", n.bootstrap))
  }
  string <- outer(s, b, paste, sep = " ") %>% 
    as.vector() %>% gsub(" +", " ", .)
  return(string)
}

phyml.sum <- function(dir) {
  files <- Sys.glob(paste0(dir, "/*")) %>% basename()
  # backup:
  system(paste0("cd ", dir, " && mkdir -p all && mv * all/"))
  # summarize by project
  trash <- files %>% split(., f = sub("_phyml_.+$", "", .)) %>% 
    lapply(function(rootname) {
      subdir <- paste0(dir, "/", rootname, "/")
      system(paste0("mkdir -p ", subdir))
      system(paste0("cd ", subdir, " && ln -s ../all/", rootname, "* ./"))
      return()
    })
  trash <- split(files, f = stringr::str_extract(files, "_phyml_.+$") %>% 
                   `[<-`(., which(is.na(.)), "input")) %>% 
    lapply(function(x) {
      if (grepl("_phyml_", x[1])) {
        subdir <- stringr::str_extract(x[1], "phyml_.+$")
      } else {
        subdir <- "input"
      }
      system(paste0("mkdir -p ", subdir, "/"))
      lapply(x, function(y) {
        system(paste0(" cd ", subdir, " && ln -s ../all/", y, " ./"))
      })
    })
  return()
}
# phyml.parse.core <- function(file) {
#   
#   browser()
#   print("m")
#   print("m")
# }
phyml.slide <- function(abao, region.gr = NULL, cons.name, max.gap.frac = 0.3,
                        window = NULL, step, threads = 1,
                        regionIDs.include = NULL, regionIDs.exclude = NULL,
                        out.dir, root.name = NULL) {
  # region.gr: locus + regionID
  if (is.character(abao)) {
    abao <- readRDS(abao)
  }
  if (is.null(abao@meta.data$consensus)) {
    abao <- aba.add.consensus(abao)
  }
  if (is.null(root.name)) {
    root.name <- abao@meta.data$fa.root.name %>% sub("_$", "", .)
  }
  # note slidingwindow returns a grangeList object
  if (is.null(region.gr)) {
    region.gr <- abao@ori.gr[1]
    all.flag <- T
  } else {
    all.flag <- F
  }
  if (is.null(window)) {
    slide.df <- region.gr %>% as.data.frame()
  } else {
    slide.gr <- region.gr %>% slidingWindows(width = window, step = step) %>% .[[1]]
    slide.df <- slide.gr %>% utilsFanc::gr2df(F, F)
  }
  
  
  if (is.null(threads)) {
    threads <- nrow(slide.df)
  }
  utilsFanc::safelapply(1:nrow(slide.df), function(i) {
    smartcut <- slide.df[i, ]
    id <- paste0("R", i, "_", utilsFanc::gr.get.loci(smartcut)) %>% sub(":", "_", .)
    if (all.flag) {
      id <- "all"
    }
    smartcut$int.name <- paste0("R", i)
    smartcut$regionID <- region.gr$regionID[1]
    smartcut$buffer.left <- 0
    smartcut$buffer.right <- 0
    fa <- paste0(out.dir, "/", root.name, "_", id, ".fa")
    phy <- sub(".fa$", ".phy", fa)
    aba.write.consensus.with.aln.2(abao = abao, cons.name = cons.name, shrink = T, 
                                   regionIDs.include = regionIDs.include, regionIDs.exclude = regionIDs.exclude,
                                   smart.cut.df = smartcut, max.gap.frac = max.gap.frac, omit.cons = T,
                                   out.file = fa)
    
    fasta2phylip(fa)
    phyml.GTR.GIF(in.phy = phy)
  }, threads = threads)
}




fedrigo.prep <- function(abao.Rds, smart.cut.df, regionIDs.include = NULL, regionIDs.exclude = NULL, 
                         out.dir, phyml = F, fedrigo.dir = NULL) {
  root.name <- basename(out.dir)
  trash <- aba.subset(abao.Rds,
                      smart.cut.df = smart.cut.df,
                      regionIDs.exclude = regionIDs.exclude, regionIDs = regionIDs.include,
                      new.work.dir = out.dir)
  aligned.fa <- paste0(out.dir, "/", root.name, "_aligned.fa")
  fa.rm.dot(aligned.fa)
  fa.aln.remove.gap(aligned.fa)
  noGap.fa <- paste0(out.dir, "/", root.name, "_aligned_noGap.fa")
  if (phyml) {
    noGap.phy <- paste0(out.dir, "/", root.name, "_aligned_noGap.phy")
    phyml.GTR.GIF(noGap.phy)
  }
  if (!is.null(fedrigo.dir)) {
    system(paste0("ln -s ", normalizePath(noGap.fa), " ", fedrigo.dir, "/"))
  }
  return()
  
}

# fedrigo.prep.tree <- function(work.dir, fedrigo.dir = NULL) {
#   root.name <- basename(work.dir)
#   tree.in <- paste0(work.dir, "/phyml/", root.name, "_aligned_noGap.phy_phyml_tree.txt")
#   tree.in <- paste0(work.dir, "/phyml/", root.name, "_aligned_noGap.phy_phyml_tree_labeled_noLength.txt")
#   
# }
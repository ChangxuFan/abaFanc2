SAMTOOLS = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
BWA = "~/software/bwa/bwa_0.7.17-r1198/bwa"
BWA.TARGET <- "~/software/bwa/bwa_0.7.16a-r1181/bwa-0.7.16/bwa"
BOWTIE2 <- "/opt/apps/bowtie2/2.3.4.1/bowtie2"
BOWTIE2.BUILD <- "/opt/apps/bowtie2/2.3.4.1/bowtie2-build"

PAFTOOLS.JS <-  "/bar/cfan/software/minimap2/misc/paftools.js"


MAFFT.AUTO <- " --auto --reorder --treeout --distout"
MAFFT.DEFAULT <- "/opt/apps/mafft/7.427/mafft"

PRANK <- "/bar/cfan/software/prank/prank/bin/prank"

BEDTOOLS <- "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools"

BEDTOOLS.DIR <- "/bar/cfan/anaconda2/envs/jupyter/bin/"

PHYML = "/bar/cfan/software/PhyML-3.1/phyml"

PHYML.params <- list()
PHYML.params$ez.SPR <- " -s SPR --rand_start --n_rand_starts 5 -b -4 "
PHYML.params$ez.GTR <- "-m GTR -c 4 -a e -v e "

PHYML.params$models <- phyml.build.models()
PHYML.params$tree_search <- phyml.build.trees()
PHYML.params$tree_search_bootstrap <- phyml.build.trees(100)

LIFTOVER_UCSC <- "/opt/apps/kentUCSC/334/liftOver"

TABIX <- "/opt/apps/htslib/1.3.1/tabix"

FINDMOTIFS.PL <- "/bar/cfan/software/homer/bin/findMotifs.pl"
HOMER.BIN <- "/bar/cfan/software/homer/bin"

FIMO <- "/opt/apps/meme/5.0.3/bin/fimo"
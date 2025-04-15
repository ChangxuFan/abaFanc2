hugo.map.to.latest <- function(df, hugo.tsv, add.cat.name = T) {
  # map df$gene to the latest hugo gene name
  # For a gene name, it will stay if it's found in the Approved.symbol
  # if not, we test if it can be found in Previous.symbols
  # if so, we replace that with Aproved.symbol
  # # if 2 Approved.symbol found for the Previous.symbol, we duplicate the row
  # else we discard the row
  utilsFanc::check.intersect("gene", "gene", colnames(df), "colnames(df)")
  
  hugo <- readr::read_tsv(hugo.tsv)
  colnames(hugo) <- make.names(colnames(hugo))
  hugo <- as.data.frame(hugo)
  
  hugo <- hugo %>% filter(Status == "Approved")
  
  if (sum(duplicated(hugo$Approved.symbol)) > 0)
    stop("sum(duplicated(hugo$Approved.symbol)) > 0")
  
  hugo2 <- tidyr::separate_rows(hugo, Previous.symbols, sep = ", ")
  hugo2 <- hugo2 %>% filter(!is.na(Previous.symbols))
  
  if (sum(is.na(hugo2$Approved.symbol)) > 0) {
    stop("sum(is.na(hugo2$Approved.symbol)) > 0")
  }
  
  bApproved <- df$gene %in% hugo$Approved.symbol
  df.app <- df[bApproved,]
  df.nApp <- df[-bApproved,]
  
  if (add.cat.name) {
    df.app$cat <- "latest"
  }
  
  bPre <- df.nApp$gene %in% hugo2$Previous.symbols
  
  if (sum(bPre) < 1) {
    return(df.app)
  }
  
  df.pre <- df.nApp[bPre,]
  
  hugo2 <- hugo2[, c("Approved.symbol", "Previous.symbols")]
  colnames(hugo2) <- c("gene.new", "gene")
  
  df.pre <- dplyr::left_join(df.pre, hugo2, by = "gene", 
                             relationship = "many-to-many")
  if (add.cat.name) {
    df.pre$cat <- paste0("P:", df.pre$gene)
  }
  
  df.pre$gene <- df.pre$gene.new
  df.pre$gene.new <- NULL
  
  df <- rbind(df.app, df.pre)
  return(df)
}

hugo.is.pseudogene <- function(x, hugo.tsv) {
  # x <- c("CRTAC1", "HSP90B2P")
  
  hugo <- readr::read_tsv(hugo.tsv)
  colnames(hugo) <- make.names(colnames(hugo))
  hugo <- as.data.frame(hugo)
  
  hugo <- hugo %>% filter(Status == "Approved")
  
  if (sum(duplicated(hugo$Approved.symbol)) > 0)
    stop("sum(duplicated(hugo$Approved.symbol)) > 0")
  
  hugo <- hugo[, c("Approved.symbol", "Approved.name")]
  colnames(hugo) <- c("symbol", "name")
  
  utilsFanc::check.intersect(x, "x", hugo$symbol, "hugo$symbol")
  
  df <- data.frame(symbol = x) %>% dplyr::left_join(hugo, by = "symbol")
  df$is.pseudo <- grepl("pseudogene", df$name)
  return(df$is.pseudo)
}

# suppressMessages(library(biomaRt))
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# gene_list <- "FAM49A"
# 
# t <- getBM(
#   attributes = c("hgnc_symbol", "external_synonym"),
#   filters = "external_synonym",
#   values = "CD8",
#   mart = ensembl
# )
# 
# mat <- read.csv(csvs[1])
# old.gn <- rownames(mat)
# 
# new.gn <- getBM(
#   attributes = c("hgnc_symbol", "external_synonym"),
#   filters = "external_synonym",
#   values = old.gn,
#   mart = ensembl
# )
# 
# new.gn %>% head()
# listFilters(ensembl) %>% View()
# 
# new.gn <- getBM(
#   attributes = c("hgnc_symbol", "external_synonym"),
#   filters = "external_synonym",
#   values = old.gn,
#   mart = ensembl
# ) 
# 
# head(new.gn)
# sum(duplicated(new.gn$external_synonym)) # 104
# dups <- new.gn$external_synonym[duplicated(new.gn$external_synonym)] %>% unique()
# 
# new.gn %>% filter(external_synonym %in% dups) %>% View()
# 
# # we use 2 step methods: if we can find the gene name in the latest version of HUGO, we 
# # will simply use it.
# 
# new.gn.pass1 <- getBM(
#   attributes = c("hgnc_symbol"),
#   filters = "hgnc_symbol",
#   values = old.gn,
#   mart = ensembl
# ) 
# 
# any(duplicated(new.gn.pass1$hgnc_symbol)) # FALSE
# 
# not.found <- old.gn %>% .[!. %in% new.gn.pass1$hgnc_symbol]
# length(not.found) # 13505
# length(old.gn) # 32738
# 
# # Half of them were not found!
# 
# new.gn.pass2 <- getBM(
#   attributes = c("hgnc_symbol", "external_synonym"),
#   filters = "external_synonym",
#   values = not.found,
#   mart = ensembl
# ) 
# 
# not.found.2 <- not.found %>% .[!.%in% new.gn.pass2$external_synonym]
# 
# length(not.found.2) # 12058
# 
# # The vast majority are just "not found"!!
# 
# not.found.2 %>% head()
# 
# hugo.map.to.latest <- function() {
#   
# }
# 



# library(org.Hs.eg.db)
# gene_list <- c("FAM49A", "FAM49B")
# # Convert to HGNC symbols
# mapped_genes <- mapIds(org.Hs.eg.db, keys = gene_list, column = "SYMBOL", keytype = "SYMBOL")
# print(mapped_genes)
# 
# ?mapIds
# columns(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)
# 
# library(httr)
# library(jsonlite)
# 
# url <- paste0("https://rest.genenames.org/search/", "FAM49A")
# 
# # Send a GET request to the API
# response <- GET(url, add_headers(Accept = "application/json"))
# 
# data <- fromJSON(content(response, "text"))
# 
# data$response$numFound

ensembl.to.symbol.gtf <- function(x, gtf, rm.version = T, must.all.map = T) {
  if (is.character(gtf))
    gr <- rtracklayer::import(gtf)
  else
    gr <- gtf
  
  df <- mcols(gr)[, c("gene_id", "gene_name")] %>% as.data.frame() %>% unique()
  colnames(df) <- c("GENEID", "SYMBOL")
  
  if (rm.version)
    df$GENEID <- df$GENEID %>% sub("\\..+$", "", .)
  
  df <- unique(df)
  if (sum(duplicated(df$GENEID)) > 0) {
    stop("sum(duplicated(df$GENEID)) > 0")
  }
  
  if (must.all.map) {
    utilsFanc::check.intersect(x, "x", df$GENEID, "df$GENEID")
  }
  
  df2 <- data.frame(GENEID = x)
  df2 <- dplyr::left_join(df2, df, by = "GENEID")
  return(df2$SYMBOL)
}

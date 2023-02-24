refbed.jsongen <- function(file, track.name = NULL, color.map.df, out.json = NULL) {
  # color.map.df: mapping the categories of column 9 in a refbed file to colors.
  # eg: data.frame(V9 = c("acti", "inhi"), color = c("black", "white"))
  # colors can also be specified by hex
  utilsFanc::check.intersect(c("V9", "color"), "required columns",
                             colnames(color.map.df), "colnames(color.map.df)")

  url <- utilsFanc::bash2ftp(file)
  if (is.null(track.name)) {
    track.name = basename(file)
  }
  if (is.null(out.json)) {
    out.json <- paste0(tools::file_path_sans_ext(sub("\\.gz$","",file)), ".json")
  }

  color.list <- as.list(color.map.df$color)
  names(color.list) <- color.map.df$V9

  jl <- list(
    name = track.name,
    type = "refbed",
    url = url,
    options = list(categoryColors = color.list)
  )
  jl <- list(jl)
  json <- jl %>% jsonlite::toJSON(auto_unbox = T) %>% jsonlite::prettify()
  dir.create(dirname(out.json), showWarnings = F, recursive = T)
  write(json, out.json)
  cat(utilsFanc::bash2ftp(out.json))
  cat("\n")
  return(json)
}

igv.hubgen <- function(files.vec, urls.vec = NULL, names.vec=NULL, xml.file=NULL, no.writing = F,
                       default.registry = "/bar/cfan/software/igv/server/default_registry.txt",
                       default.xml = "/bar/cfan/software/igv/server/default_xml.txt",
                       global.name = NULL, registry.file = NULL, override.xml, override.registry) {
  if (is.null(urls.vec))
    urls.vec <- files.vec %>% utilsFanc::bash2ftp()
  if (is.null(names.vec)) 
    names.vec <- basename(urls.vec)
  if (is.null(global.name))
    global.name <- "fanc"
  node <- XML::newXMLNode("Global", attrs = list(name = global.name))
  mapply(function(url, name) {
    XML::newXMLNode("Resource", attrs = list(name = name, path = url), parent = node)
  }, urls.vec, names.vec, SIMPLIFY = F)
  xml <- XML::saveXML(node)
  
  if (no.writing == T) 
    return(cat(xml))
  
  if (is.null(xml.file)) {
    xml.file <- default.xml
    override.xml <- T
  }
    
  if (is.null(registry.file)) {
    registry.file <- default.registry
    override.registry <- T
  }
  
  
  if (override.xml == T) {
    try(system(paste0("mv ", default.xml, " ", default.xml, ".bk")))
    write(xml, file = xml.file)
    if (xml.file != default.xml)
      try(system(paste0("cp ", xml.file, " ", default.xml)))
  }
  
  if (override.registry == T) {
    try(system(paste0("mv ", default.registry, " ", default.registry, ".bk")))
    write(xml.file  %>% utilsFanc::bash2ftp(), file = registry.file)
    if (registry.file != default.registry)
      try(system(paste0("cp ", registry.file, " ", default.registry)))
  }
  
  return(cat(xml))
}


dir.2.igv <- function(dir, pattern = NULL, recursive, ...) {
  files <- list.files(path = dir, pattern = NULL, full.names = F, recursive = recursive,
                      no.. = T, all.files = F, ignore.case = F, include.dirs = T)
  if (!is.null(pattern))
    files <- files[grepl(pattern, files)]
  files <- paste0(dir, "/", files) %>% normalizePath(mustWork = T)
  igv.hubgen(files, ...) %>% 
  return()
}



### xref functions


process_xref <- function(xfile){
  # foreign <- read.table(xfile,
  #                       header = TRUE,
  #                       sep = "\t",
  #                       stringsAsFactors = FALSE)
  ## Modified by Sergio Romera: adding quote and fill to read.table to deal with missing values in dme ref
  foreign <- read.table(xfile,
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        quote="",
                        fill=FALSE)
  this_entrezs <- unique(foreign[,1])
  foreign <- foreign[which(foreign[,2]!=""),]
  foreign <- unique(foreign)
  xref <- by(foreign, foreign[,2], function(x) x[,1])
  s <- summary(sapply(xref, length))
  tt <- table(sapply(xref, length))
  return(list(xref = xref,
              entrezs = this_entrezs,
              n = nrow(foreign),
              name = xfile,
              stats = cbind(do.call("cbind",as.list(s)),
                            n = nrow(foreign)),
              freqs = tt))
}

process_species <- function(xref_folder){
  xref_files <- list.files(xref_folder, full.names = TRUE)
  xref_rdata <- paste0(xref_folder, "/xref_", xref_folder, ".RData")
  xref_files <- xref_files[grep(xref_rdata, xref_files, invert = TRUE)]
  xref <- list()
  all_entrezs <- c()
  stats <- c()
  tts <- list()
  for(xfile in xref_files){
    cat("Processing", xfile, "\n")
    this_xref <- process_xref(xfile)
    print(this_xref$stats)
    xref <- c(xref, this_xref$xref)
    all_entrezs <- c(all_entrezs, this_xref$entrezs)
    stats <- rbind(stats, this_xref$stats)
    tts[[xfile]] <- this_xref$freqs
  }
  all_entrezs <- unique(all_entrezs)
  entrezs_list <- as.list(all_entrezs)
  names(entrezs_list) <- all_entrezs
  xref <- c(xref, entrezs_list)
  rownames(stats) <- xref_files
  attr(xref, "stats") <- stats
  attr(xref, "freqs") <- tts
  save(xref, file = paste0(xref_folder, "/xref_", xref_folder, ".RData"))
}

### process species

#setwd("RDatas_constructor/xref")
#### Modified by Sergio Romera: Adding new path to xref refs subdirectory
setwd("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/xref/refs/")
#species <- c("hsa", "mmu", "rno")
#### Modified by Sergio Romera: Adding dme as new species
species <- c("hsa", "mmu", "rno", "dme")

for(spe in species){
  cat(">>> Processing", spe, "\n")
  process_species(spe)
  cat("\n\n")
}


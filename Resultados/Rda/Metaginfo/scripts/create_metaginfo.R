## marta.hidalgo@outlook.es
## Create metaginfo objects 

rm(list = ls())

library(KEGGgraph)
library(igraph)
library(graph)
# library(hipathia)
## Modified by Sergio Romera: loading local hipathia package
library(devtools)
load_all("/home/sromera/Documents/TFM/scripts_hipathia/ultimate_hipathia/hipathia-master/")
load_all("/home/sromera/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpAnnot/")
# library(hpAnnot)

# path <- "~/appl/hpConstructor/RDatas_constructor/pathways/"
path <- "~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/pathways/"
setwd(path)
hipath <- paste0(path, "scripts/")
# source(paste0(hipath, "/graphs.R"))
# source(paste0(hipath, "/KEGG_net.R"))
# source(paste0(hipath, "/layout.R"))
source(paste0(hipath, "graphs.R"))
source(paste0(hipath, "KEGG_net.R"))
source(paste0(hipath, "layout.R"))
# source("~/appl/hipathia/R/utils.R")
source("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hipathia-master/R/utils.R")
source("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hipathia-master/R/load_hpAnnot.R")

# ammend.file <- paste0(hipath, "/../sif_amendments.txt")
# comp.file <- paste0(hipath, "/../compounds_list.txt")
ammend.file <- paste0(path, "sif_amendments.txt")
comp.file <- paste0(path, "compounds_list.txt")

# FUNCTION
load_compounds <- function(comp.file){
    compounds <- utils::read.table(comp.file, header = FALSE, sep = "\t", 
                                   stringsAsFactors = FALSE, 
                                   colClasses = "character", quote = "")
    colnames(compounds) <- c("ID", "names")
    compounds$simple.ID <- gsub("cpd:", "", compounds$ID)
    compounds$simple.name <- sapply(strsplit(compounds$names, split = ";"), "[[", 1) 
    return(compounds[,c(3,4)])
}

compounds <- load_compounds(comp.file)
# save(compounds, file=paste0(hipath, "/../compounds_list.RData"))

# Parameters
# species <- c("hsa", "rno", "mmu")
species <- c("hsa", "rno", "mmu", "dme")
date <- gsub("-", "_", Sys.Date())
## Modified by Sergio Romera
version <- "v2"

# species <- "dme"
# spe <- "rno"
# spe <- "hsa"
spe <- "dme"
## Modified by Sergio Romera
# path <- "~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/pathways"
for(spe in species){

    # set folders
    kgml.folder <- paste0(path, spe, "/kgml/")
    sif.folder <- paste0(path, spe, "/sif/")
    tmp.folder <- paste0(path, spe, "/temp/")
    pathway.names <- unique(gsub(".xml", "", list.files(kgml.folder, 
                                                        pattern="xml")))
    
    # Create folders
    if(!dir.exists(sif.folder))
        dir.create(sif.folder)
    if(!dir.exists(tmp.folder))
        dir.create(tmp.folder)
    
    # Load annotations
    # dbannot <- hipathia:::load_annots("uniprot", spe)
    # entrez2hgnc <- hipathia:::load_entrez_hgnc(spe)
    # Modified by Sergio Romera: it loads the functions from 
    dbannot <- load_annots("uniprot", spe)
    entrez2hgnc <- load_entrez_hgnc(spe)

    # Process KGML files
    #-------------------------------------------------
    # create name_pathways.txt file
    create.pathway.names.file(pathway.names, spe, kgml.folder, sif.folder)

    # Transform KGML to SIF files
    transform.XML.to.SIF(pathway.names, kgml.folder, sif.folder, compounds)

    # Load pathways from created SIF files
    pgs <- load.graphs(sif.folder, spe)
    save(pgs, file=paste0(tmp.folder, "/pgs_", date, ".RData"))

    # Ammend pathways
    apgs <- amend.kegg.pathways(ammend.file, pgs, spe)
    save(apgs, file=paste0(tmp.folder, "/apgs_", date, ".RData"))

    # Add final functions to the pathways
    fpgs <- add.functions.to.pathigraphs(apgs, entrez2hgnc, dbannot,
                                         maxiter = 1000)
    save(fpgs, file=paste0(tmp.folder, "/fpgs_", date, ".RData"))
    # load(file=paste0(tmp.folder, "/fpgs_", date, ".RData"))
    
    # Compute Path Normalization Values
    metaginfo <- create.metaginfo.object(fpgs, spe)
    save(metaginfo, file=paste0(tmp.folder, "/meta_graph_info_", spe,
                                "_", date, ".RData"))

}



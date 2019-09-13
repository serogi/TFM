
library(igraph)
library(devtools)
# library(hipathia)
# library(hpAnnot)
## Modified by Sergio Romera
load_all("/home/sromera/Documents/TFM/scripts_hipathia/ultimate_hipathia/hipathia-master/")
load_all("/home/sromera/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpAnnot/")

species <- c("hsa", "mmu", "rno", "dme")
spe <- "dme"
# rda_path <- "RDatas_constructor/RDatas/"
## Modified by Sergio Romera
rda_path <- "~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor"
edpath <- "/home/sromera/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpAnnot/extdata/"
date <- gsub("-", "_", Sys.Date())
version <- "v2"

for(spe in species){

    print(spe)
    # mgi <- load(paste0(rda_path, "/meta_graph_info_", spe, "_", version, ".rda"))
    mgi <- load(paste0(rda_path,"/pathways/XMLs/", spe, "/temp/meta_graph_info_", spe, "_", date, ".RData"))
    metaginfo <- get(mgi)
    
    # Uniprot
    # unidb <- load(paste0(rda_path,  "annot_uniprot_", spe, "_", version, ".rda"))
    ## Modified by Sergio Romera
    unidb <- load(paste0(edpath, "v2/annot_uniprot_", spe, "_", version, ".rda"))
    uni_bp_annot <- get(unidb)
    # annofuns <- hipathia:::annotate_paths(metaginfo, uni_bp_annot)
    ## Modified by Sergio Romera
    annofuns <- annotate_paths(metaginfo, uni_bp_annot)
    # save(annofuns, file = paste0("RDatas_constructor/annofuns/annofuns_uniprot_",
    #                              spe, ".RData"))
    ## Modified by Sergio Romera
    save(annofuns, file = paste0("annofuns/annofuns_uniprot_",
                                 spe, ".RData"))
    

    # GO
    # godb <- load(paste0(rda_path, "annot_GO_", spe, "_", version, ".rda"))
    ## Modified by Sergio Romera
    godb <- load(paste0(edpath, "/v2/annot_GO_", spe, "_", version, ".rda"))
    go_bp_annot <- get(godb)
    # annofuns <- hipathia:::annotate_paths(metaginfo, go_bp_annot)
    ## Modified by Sergio Romera
    annofuns <- annotate_paths(metaginfo, go_bp_annot)
    # save(annofuns, file = paste0("RDatas_constructor/annofuns/annofuns_GO_",
    #                              spe, ".RData"))
    ## Modified by Sergio Romera
    save(annofuns, file = paste0("annofuns/annofuns_GO_",
                                 spe, ".RData"))


}


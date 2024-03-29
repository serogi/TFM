

# species <- "hsa"
# species <- "rno"
# species <- "mmu"
# species <- "dme"

# all_species <- c("hsa", "mmu", "rno")
## Modified by Sergio R: added dme to species
all_species <- c("hsa", "mmu", "rno", "dme")
# dir.create("RDatas_constructor/annotations/annotations/")
## Modified by Sergio R: added path to work with local hipathia files
dir.create("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/annotations/")

for(species in all_species){

    print(species)

    # ann_path <- "RDatas_constructor/annotations/annotations/"
    # raw_path <- "RDatas_constructor/annotations/raw_data/"
    # ann_spe_path <- paste0("RDatas_constructor/annotations/annotations/", species, "/")
    # raw_spe_path <- paste0("RDatas_constructor/annotations/raw_data/", species, "/")
    ann_path <- "~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/annotations/"
    raw_path <- "~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/raw_data/"
    ann_spe_path <- paste0("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/annotations/", species, "/")
    raw_spe_path <- paste0("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/raw_data/", species, "/")
    if(!dir.exists(paste0("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/annotations/", species))){
        dir.create(paste0("~/Documents/TFM/scripts_hipathia/ultimate_hipathia/hpConstructor-master/RDatas_constructor/annotations/annotations/", species))
    }

    
    ## hgnc to entrez
    rh_file <- paste0(raw_spe_path, "/raw_entrez_hgnc_", species, ".txt")
    raw_hgnc <- read.delim(rh_file,
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE)
    
    # if(species == "hsa"){
    #     raw_hgnc <- raw_hgnc[,c("entrez","hgnc")]
    #     colnames(raw_hgnc) <- c("EntrezGene.ID", "Associated.Gene.Name")
    # }
    # # Modified by Sergio Romera: changing colnames in rno and dme
    if(species == "hsa"){
      raw_hgnc <- raw_hgnc[,c("entrez","hgnc")]
      colnames(raw_hgnc) <- c("EntrezGene.ID", "Associated.Gene.Name")
    } else if(species == "dme"){
      raw_hgnc <- raw_hgnc[,c("EntrezGene.ID", "Gene.name")]
      colnames(raw_hgnc) <- c("EntrezGene.ID", "Associated.Gene.Name")
    } else if(species == "rno"){
      raw_hgnc <- raw_hgnc[,c("NCBI.gene.ID", "Gene.name")]
      colnames(raw_hgnc) <- c("EntrezGene.ID", "Associated.Gene.Name")
    }
    clean_entrez_hgnc <- raw_hgnc[!is.na(raw_hgnc$EntrezGene.ID),]
    agn <- clean_entrez_hgnc$Associated.Gene.Name
    selcols <- c("EntrezGene.ID", "Associated.Gene.Name")
    clean_entrez_hgnc <- clean_entrez_hgnc[!is.na(agn),selcols]
    eh_file <- paste0(ann_spe_path, "/entrez_hgnc_", species, ".annot")
    write.table(clean_entrez_hgnc,
                file = eh_file,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE)
    
    ## uniprot keywords
    ru_file <- paste0(raw_spe_path, "uniprot_", species, "_keywords_noprime.annots")
    raw_uniprot <- read.delim(ru_file,
                              sep="\t",
                              header = TRUE,
                              stringsAsFactors = FALSE)

    # create_gene_annotations <- function(x){
    #     genes <- unlist(strsplit(as.character(x['Gene.names'])," "))
    #     x['Keywords'] <- gsub("; ", ";", x['Keywords'])
    #     keywords <- unlist(strsplit(as.character(x['Keywords']),";"))
    #     annots <- cbind(rep(genes, each = length(keywords)),
    #                     rep(keywords, times = length(genes)))
    #     return(annots)
    # }
    # uniprot_keywords <- do.call("rbind", apply(raw_uniprot, 1,
    #                                            create_gene_annotations))
    ## Modified by Sergio Romera: changing the create_gene_annotations step to split dme keywords
    ## by ";" instead of "; "
    if(species == "dme"){
      
      create_gene_annotations_dme <- function(x){
        genes <- unlist(strsplit(as.character(x['Gene.names'])," "))
        keywords <- unlist(strsplit(as.character(x['Keywords']),";")) # Well, hello there
        annots <- cbind(rep(genes, each = length(keywords)),
                        rep(keywords, times = length(genes)))
        return(annots)
      }
      uniprot_keywords <- do.call("rbind", apply(raw_uniprot, 1,
                                                 create_gene_annotations_dme))
    }else{
      create_gene_annotations <- function(x){
        genes <- unlist(strsplit(as.character(x['Gene.names'])," "))
        keywords <- unlist(strsplit(as.character(x['Keywords']),"; "))
        annots <- cbind(rep(genes, each = length(keywords)),
                        rep(keywords, times = length(genes)))
        return(annots)
      }
      uniprot_keywords <- do.call("rbind", apply(raw_uniprot, 1,
                                                 create_gene_annotations))
    }

    require("XML")
    xkey <- xmlParse(paste0(raw_path, "/keywords.rdf"))
    xmltop = xmlRoot(xkey)
    xmlList <- xmlToList(xmltop)
    parse_node <- function(node){
        fields <- names(node)
        name <- node$prefLabel
        id <- basename(as.character(node[['.attrs']]))
        parent <- basename(as.character(node[which(fields=="subClassOf")]))
        if(is.null(name)){
            return(NULL)
        } else {
            id <- paste("n", id, sep = "_")
            if(length(parent) == 0){
                parent <- "root"
            } else {
                parent <- paste("n", parent, sep = "_", collapse = ",")
            }
            c(name = name, id = id, parent = parent)
        }
    }
    keyonto <- as.data.frame(do.call("rbind",lapply(xmlList, parse_node)),
                             stringsAsFactors = FALSE)
    rownames(keyonto) <- keyonto$id

    get_super <- function(x){
        if(keyonto[x,"parent"] == "root"){
            return(x)
        } else {
            parents <- unlist(strsplit(keyonto[x,"parent"], ","))
            unique(unlist(sapply(parents, get_super)))
        }
    }
    supers <- sapply(keyonto$id, get_super)
    keyonto$super <- sapply(supers, function(x) paste(x, collapse=","))
    keyonto$super.name <- sapply(supers, function(x) paste(keyonto[x,"name"],
                                                           collapse=","))

    uniprot_keywords <- cbind(uniprot_keywords,
                              as.character(
                                  keyonto$super.name[match(uniprot_keywords[,2],
                                                           keyonto$name)]))
    uk_file <- paste0(ann_spe_path, "/uniprot_keywords_", species, ".annot")
    write.table(uniprot_keywords,
                file = uk_file,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE)

    main_cat <- keyonto[which(keyonto$parent=="root"),"name"]
    for(mc in main_cat){
        fn <- paste0(ann_spe_path, "/uniprot_keywords_", species, "__",
                     gsub(" ", "_", tolower(mc)), ".annot")
        write.table(uniprot_keywords[grep(mc, uniprot_keywords[,3]),],
                    file = fn,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = "\t",
                    quote = FALSE)
    }





    ## gene ontology
    library(limma)
    library(igraph)

    obo.parser <- function(obo_file){

        terms <- list()
        id <- "General"
        terms[[id]] <- list()
        terms[[id]]$id <- id

        con <- file(obo_file)
        open(con)
        while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
            fields <- unlist(strsplit(line, ": "))
            if(length(fields) > 1){
                key <- fields[1]
                value <- fields[2]
                if(key=="id"){
                    id <- trimWhiteSpace(value)
                    terms[[id]] <- list()
                    terms[[id]]$id <- id
                }
                else {
                    if(is.null(terms[[id]][[key]])){
                        terms[[id]][[key]] <- value
                    } else {
                        terms[[id]][[key]] <- c(terms[[id]][[key]], value)
                    }
                }
            }
        }
        close(con)
        general <- terms$General
        terms$General <- NULL
        attr(terms,"General") <- general
        return(terms)
    }

    # #gos <- obo.parser("test.obo")
    # gos <- obo.parser("RDatas_constructor/annotations/go-basic.obo")
    # save(gos,file="RDatas_constructor/annotations/go-basic.rdata")
    load(paste0(raw_path, "/go-basic.rdata"))

    go_namespace <- sapply(gos, "[[", "namespace")

    get.go.frame <- function(gosel){
        go_frame <- do.call("rbind", lapply(gosel, function(x) {
            if(is.null(x$is_a)){
                parents <- ""
            } else {
                parents <- paste(sapply(strsplit(x$is_a," ! "), "[[", 1),
                                 collapse = ",")
            }
            return(c(id = x$id,
                     name = x$name,
                     namespace = x$namespace,
                     parents = parents))
        }))
        return(as.data.frame(go_frame, stringsAsFactors = FALSE))
    }

    get.go.net <- function(gosel){
        go_edges <- do.call("rbind", sapply(gosel, function(x) {
            if(is.null(x$is_a)){
                return(NULL)
            } else {
                parents <- sapply(strsplit(x$is_a, " ! "), "[[", 1)
                return(cbind(parents, rep(x$id, length(parents))))
            }
        }))
        go_net <- graph.data.frame(go_edges, directed = TRUE)
        return(go_net)
    }

    gos_bp <- gos[go_namespace == "biological_process"]
    go_bp_frame <- get.go.frame(gos_bp)
    go_bp_net <- get.go.net(gos_bp)
    dis <- shortest.paths(go_bp_net, v = "GO:0008150")[1,,drop = TRUE]+1
    go_bp_frame$level <- dis[match(rownames(go_bp_frame),names(dis))]


    raw_go_file <- paste0(raw_spe_path, "/gene_GO_association_", species, ".txt")
    raw_go_annots <- read.delim(raw_go_file,
                                sep = "\t",
                                header = FALSE,
                                stringsAsFactors = FALSE,
                                comment.char = "!")
    go_bp_annots <- raw_go_annots[,c(3,5,7,9)]
    colnames(go_bp_annots) <- c("gene", "term", "evidence", "namespace")
    go_bp_annots <- go_bp_annots[ go_bp_annots$namespace == "P", ]

    gbf_file <- paste0(ann_spe_path, "go_bp_frame_", species, ".txt")
    write.table(go_bp_frame,
                file = gbf_file,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)
    gba_file <- paste0(ann_spe_path, "go_bp_", species, "_annots.txt")
    write.table(go_bp_annots[,1:3],
                file = gba_file,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)

    save(go_bp_net, file = paste0(ann_path, "/go_bp_net.RData"))
    save(go_bp_frame, file = paste0(ann_path, "/go_bp_frame.RData"))

    minigo <- go_bp_annots[go_bp_annots$evidence == "EXP" |
                               go_bp_annots$evidence == "IDA", ]
    namegos <- go_bp_frame[minigo[,2],2]
    goname <- as.data.frame(cbind(gene = minigo$gene, "function" = minigo$term, term = namegos))
    goname_file <- paste0(ann_spe_path, "go_bp_", species, ".annot")
    write.table(goname,
                file = goname_file,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)


}






# Hasta aquí lo que usamos de normal

# ##### go slim
# goslim <- obo.parser("goslim_generic.obo")
# goslim_frame <- get.go.frame(goslim)
# goslim_bp_frame <- goslim_frame[ goslim_frame$namespace=="biological_process",]
# gos_bp_frame <- read.delim("go_bp_frame.txt",header=T,sep="\t",stringsAsFactors=F)
# annots <- read.delim("go_bp_human_annots.txt",header=T,sep="\t",stringsAsFactors=F)
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# annots_slim <- annots[ annots$term %in% goslim_bp_frame$id,]
# write.table(annots_slim,file="go_slim_bp_human_annots.txt",sep="\t",col.names=T,row.names=F,quote=F)
#
#
# # others
# gonot <- obo.parser("gocheck_do_not_annotate.obo")
# gonot_frame <- get.go.frame(gonot)
# gonot_bp_frame <- gonot_frame[ gonot_frame$namespace=="biological_process",]
# gonotman <- obo.parser("gocheck_do_not_manually_annotate.obo")
# gonotman_frame <- get.go.frame(gonotman)
# gonotman_bp_frame <- gonotman_frame[ gonotman_frame$namespace=="biological_process",]
#
# # test
# # load files
#
# gene <- "BAD"
# gene <- "LIPE"
# gene <- "PYGL"
# gene <- "FBP1"
# gene <- "RPS6"
# gene <- "PYGB"
#
# gene <- "CASP3"
#
# gene <- "CDK4"
# gene <- "CDKN1A"
# gene <- "SGK2"
# gene <- "PKN1"
# gene <- "PRKCA"
# gene <- "BRCA1"
# gene <- "FBP1"
# gene <- "FASN"
# gene <- "GYS2"
# gene <- "EIF4B"
#
#   # go
# invalid_evidences <- "IEA"
# annots <- go_bp_annots[-which(go_bp_annots$evidence %in% invalid_evidences),]
# valid_evidences <- c("EXP","IDA","IPI","IMP","IGI","IEP","TAS","ISS")
# annots <- go_bp_annots[which(go_bp_annots$evidence %in% valid_evidences),]
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# #annots <- annots[-grep("signaling pathway",annots$name),]
# annots <- unique(annots)
#
# gene_annots <- annots[ annots$gene==gene, ]
# dups <- which(duplicated(gene_annots$term)==T)
# if(length(dups)>0) gene_annots <- gene_annots[-dups,]
# gene_annots <- gene_annots[order(gene_annots$level),]
# gene_annots
# #gene_net <- induced.subgraph(go_bp_net,v = gene_annots$term)
#
#
#   #goslim
# annots <- go_bp_annots
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# annots_slim <- annots[ annots$term %in% goslim_bp_frame$id,]
# gene_annots <- annots_slim[ annots_slim$gene==gene, ]
# gene_annots
#
# #signaling
# descriptions <- sapply(gos_bp,"[[","def")
# gos_bp_frame_sig <- gos_bp_frame[grep("signaling|pathway|apoptotic",
# descriptions),]
# annots <- go_bp_annots
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# annots_sig <- annots[ annots$term %in% gos_bp_frame_sig$id,]
# valid_evidences <- c("EXP","IDA","IPI","IMP","IGI","IEP","TAS","ISS")
# annots_sig <- annots_sig[which(annots_sig$evidence %in% valid_evidences),]
# gene_annots <- annots_sig[ annots_sig$gene==gene, ]
# gene_annots
#
# # others
# annots <- go_bp_annots
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# annots_gonot <- annots[ annots$term %in% gonot_bp_frame$id,]
# gene_annots <- annots_gonot[ annots_gonot$gene==gene, ]
# gene_annots
#
# annots <- go_bp_annots
# annots$name <- gos_bp_frame$name[match(annots$term,gos_bp_frame$id)]
# annots$level <- gos_bp_frame$level[match(annots$term,gos_bp_frame$id)]
# annots_gonotman <- annots[ annots$term %in% gonotman_bp_frame$id,]
# gene_annots <- annots_gonotman[ annots_gonotman$gene==gene, ]
# gene_annots
#
#
# dbannot <- annots[,c("gene","name")]
# pathfuncs <- lapply(pathigraph$subgraphs,get.path.functions,dbannot=dbannot,
# entrez2hgnc=entrez2hgnc)
#
#
# ############################### HUMAN PHENOTYPE ONTOLOGY
#
#
# hpo <- obo.parser("hp.obo")
# fields <- c("id","name","def")
# hpo_frame <- do.call("rbind",lapply(hpo,"[",fields))
# hpo_frame[which(sapply(hpo_frame[,"def"],is.null)==T),"def"] <- ""
# hpo_frame[,"def"] <- gsub("\"","",hpo_frame[,"def"])
# colnames(hpo_frame) <- fields
# write.table(hpo_frame,file="hpo_frame.txt",row.names=F,col.names=T,sep="\t",
# quote=F)
#
# hpo_annot <- read.delim("ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt",
# sep="\t",header=T,stringsAsFactors=F)
# write.table(hpo_annot[,c(2,4)],file="hpo_annot.txt",row.names=F,col.names=F,
# sep="\t",quote=F)
#
#
#
#
#
#






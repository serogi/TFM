

## GEOexpressionmatrix
## Obtiene una matriz de expresión de las muestras del experimento de GEO, con los nombres de las filas
## correspondientes a los identificadores de genes proporcionados por el usuario.

GEOexpressionmatrix <- function(geo_entry,geo_matrix_id, gene_id) {
  geo_file <- paste0(geo_matrix_id,"_series_matrix.txt.gz")
  expr_matrix <- geo_entry[[geo_file]]@assayData$exprs
  list_id <- geo_entry[[geo_file]]@featureData@data[[gene_id]]
  rownames(expr_matrix) <- list_id
  return(expr_matrix)
}


## GEOdesignmatrix
## Obtiene una matriz con los datos de las variables estudiadas en el experimento
## Es necesario añadir el número de variables encontradas

GEOdesignmatrix <- function(geo_entry, geo_matrix_id, n_characteristics){
  geo_file <- paste0(geo_matrix_id,"_series_matrix.txt.gz")
  n_samples <- length(levels(geo_entry[[geo_file]]@phenoData@data[["title"]]))
  design_matrix <- as.data.frame(matrix(ncol=n_characteristics, nrow = n_samples))
  characteristic = "characteristics_ch1"
  for(char in 1:n_characteristics){
    if(char == 1){
      design_matrix[,char] <- geo_entry[[geo_file]]@phenoData@data[[
        which(colnames(geo_entry[[geo_file]]@phenoData@data) == characteristic)]]
      colnames(design_matrix)[char] <- characteristic
    }
    else{
      characteristic2 <- paste0(characteristic,".",(char-1))
      design_matrix[,char] <- geo_entry[[geo_file]]@phenoData@data[[
        which(colnames(geo_entry[[geo_file]]@phenoData@data) == characteristic2)
        ]]
      colnames(design_matrix)[char] <- characteristic2
    }
  }
  rownames(design_matrix) <- geo_entry[[geo_file]]@phenoData@data[["geo_accession"]]
  return(design_matrix)
}

## GEO2SummarizedExperiment
## Transformamos matrices de expresión (previamente normalizada) y de diseño a un formato SummarizedExperiment

GEO2SummarizedExperiment <- function(normalized_expr_matrix, design_matrix){
  SummExp <- SummarizedExperiment(assays=SimpleList(raw=normalized_expr_matrix), colData=design_matrix)
  return(SummExp)
}


## PCA2DataFrame
## Creamos un data frame a partir de una PCA obtenida de un SummarizedExperiment (Requiere mixomics)

PCA2DataFrame <- function(summ_exp){
  pca <- mixOmics::pca(t(assay(summ_exp)), ncomp = 4)
  pc1 <- pca$x[,1]
  pc2 <- pca$x[,2]
  pc3 <- pca$x[,3]
  pc4 <- pca$x[,4]
  data_frame <- as.data.frame(summ_exp@colData@listData)
  rownames(data_frame) <- summ_exp@colData@rownames
  data_frame$PC1 <- pc1 
  data_frame$PC2 <- pc2
  data_frame$PC3 <- pc3
  data_frame$PC4 <- pc4
  return(data_frame)
}


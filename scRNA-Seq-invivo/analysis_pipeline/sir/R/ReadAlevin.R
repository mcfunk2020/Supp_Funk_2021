
#' https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/
ReadAlevinCsv <- function( base.path = NULL ){
  barcode.loc <- paste0( base.path, "alevin/quants_mat_rows.txt" )
  gene.loc <- paste0( base.path, "alevin/quants_mat_cols.txt" )
  matrix.loc <- paste0( base.path, "alevin/quants_mat.csv" )

  matrix <- as.matrix(utils::read.csv( matrix.loc, header=FALSE))
  matrix <- t(matrix[,1:ncol(matrix)-1])

  cell.names <- readLines( barcode.loc )
  gene.names <- readLines( gene.loc )

  colnames(matrix) <- cell.names
  rownames(matrix) <- gene.names
  matrix[is.na(matrix)] <- 0
  return(matrix)
}

ReadAlevinGz <- function( base.path = NULL ){
  .NotYetImplemented()
  barcode.loc <- paste0( base.path, "alevin/quants_mat_rows.txt" )
  gene.loc <- paste0( base.path, "alevin/quants_mat_cols.txt" )
  matrix.loc <- paste0( base.path, "alevin/quants_mat.csv" )

  matrix <- as.matrix(utils::read.csv( matrix.loc, header=FALSE))
  matrix <- t(matrix[,1:ncol(matrix)-1])

  cell.names <- readLines( barcode.loc )
  gene.names <- readLines( gene.loc )

  colnames(matrix) <- cell.names
  rownames(matrix) <- gene.names
  matrix[is.na(matrix)] <- 0
  return(matrix)
}

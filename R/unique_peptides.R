#' BuildPeptideProteinMatrix
#'
#' Build an adjacency matrix for peptides and proteins.
#' This function is strongly inspired by BuildAdjacencyMatrix() from the DAPAR module
#' of the Prostar package (c:\Users\adrie\Documents\LSMBO\DAPAR\R\agregation.R, line 351).
#' It uses the same mathematical approach but with sparseMatrix() instead of table()
#' to handle very large datasets (>700k precursors) efficiently.
#' where rows = peptides/precursors and columns = proteins.
#'
#' @param df Data frame containing the data (precursors)
#' @param peptide_col Name of the column containing the peptide/precursor identifier
#' @param protein_col Name of the column containing proteins (may contain ";" for shared peptides)
#'
#' @return A sparse adjacency matrix
#'
#' @export
BuildPeptideProteinMatrix <- function(df, peptide_col = "Precursor.Id", protein_col = "Protein.Group") {
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for this function.")
  }
  
  # Extract protein groups column
  PG <- df[[protein_col]]
  
  # Split protein groups and trim whitespace (inspired by DAPAR approach)
  # Creates a list where each element contains the proteins for one peptide
  PG.l <- lapply(
    strsplit(as.character(PG), "[;]+"),
    function(x) trimws(x)
  )
  
  # Build sparse adjacency matrix efficiently for large datasets
  # This approach is inspired by DAPAR's BuildAdjacencyMatrix() but adapted
  # to handle very large datasets (>700k precursors) by using sparseMatrix directly
  # instead of table() which has a 2^31 elements limitation
  
  # Get all unique proteins
  all_proteins <- unique(unlist(PG.l))
  
  # Build index vectors for sparseMatrix
  # i = row indices (peptide indices, repeated for each protein)
  # j = column indices (protein indices)
  i <- rep(seq_along(PG.l), lengths(PG.l))
  j <- match(unlist(PG.l), all_proteins)
  
  # Create sparse matrix with 1's where peptide-protein associations exist
  mat <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(PG.l), length(all_proteins)),
    dimnames = list(df[[peptide_col]], all_proteins)
  )
  
  return(mat)
}


#' DeterminePeptideUniqueness
#'
#' Determine if each peptide/precursor is unique or shared.
#' Inspired by Prostar logic: a peptide is unique if it is associated with only one protein
#' in the adjacency matrix.
#'
#' @param adjacency_matrix Adjacency matrix for peptides and proteins
#'
#' @return A named vector with TRUE for unique peptides, FALSE for shared peptides
#'
#' @export
DeterminePeptideUniqueness <- function(adjacency_matrix) {
  
  # For each row (peptide), count the number of associated proteins
  # Sum over columns: if = 1, unique peptide; if > 1, shared peptide
  nb_proteins_per_peptide <- Matrix::rowSums(adjacency_matrix)
  
  # A peptide is unique if associated with exactly 1 protein
  is_unique <- nb_proteins_per_peptide == 1
  
  return(is_unique)
}


#' AddUniquenessColumn
#'
#' Add a "Unique" column to the dataframe.
#' Main function that combines matrix construction and uniqueness determination.
#'
#' @param df Data frame containing the data
#' @param peptide_col Name of the peptide/precursor column
#' @param protein_col Name of the protein column
#'
#' @return The dataframe with an added "Unique" column ("1" for unique, "0" for shared)
#'
#' @export
AddUniquenessColumn <- function(df, 
                               peptide_col = "Precursor.Id", 
                               protein_col = "Protein.Group") {
  
  # Build adjacency matrix
  adj_matrix <- BuildPeptideProteinMatrix(df, peptide_col, protein_col)
  
  # Determine uniqueness (returns a named logical vector)
  is_unique <- DeterminePeptideUniqueness(adj_matrix)
  
  # Vectorized approach: match peptide IDs to the is_unique vector
  # Much faster than for loop for large datasets (700k+ rows)
  peptide_ids <- df[[peptide_col]]
  df$Unique <- ifelse(is_unique[peptide_ids], "1", "0")
  
  # Handle NAs (peptides not found in matrix, shouldn't happen but safety check)
  df$Unique[is.na(df$Unique)] <- "0"
  
  return(df)
}


#' SaveDataWithUnique
#'
#' Save the dataframe with Unique column as TSV file.
#'
#' @param df Dataframe to save
#' @param output_path Path to the output file
#'
#' @export
SaveDataWithUnique <- function(df, output_path) {
  
  # Save as TSV
  write.table(df, 
              file = output_path,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  message(paste0("File saved: ", output_path))
}


#' ProcessUniquePeptides
#'
#' Main function for the complete workflow.
#' This function combines all steps:
#' 1. Add Unique column
#' 2. Save debug TSV file
#' 3. Filter to keep only unique peptides if requested
#'
#' @param df Input dataframe
#' @param save_debug If TRUE, save the file with Unique column
#' @param debug_path Path to the debug file
#' @param filter_unique If TRUE, keep only rows with Unique="1"
#' @param peptide_col Name of the peptide/precursor column
#' @param protein_col Name of the protein column
#'
#' @return The dataframe (filtered or not depending on filter_unique)
#'
#' @export
ProcessUniquePeptides <- function(df,
                                 save_debug = TRUE,
                                 debug_path = file.path(getwd(), "data_with_unique_col.tsv"),
                                 filter_unique = TRUE,
                                 peptide_col = "Precursor.Id",
                                 protein_col = "Protein.Group") {
  
  # Step 1: Add Unique column
  df_with_unique <- AddUniquenessColumn(df, peptide_col, protein_col)
  
  # Step 2: Save debug file if requested
  if (save_debug) {
    SaveDataWithUnique(df_with_unique, debug_path)
  }
  
  # Step 3: Filter if requested
  if (filter_unique) {
    df_filtered <- df_with_unique[df_with_unique$Unique == "1", ]
    return(df_filtered)
  } else {
    return(df_with_unique)
  }
}


#' ApplyUniquePeptidesFilter
#'
#' Simplified compatibility function for easy integration into app.R.
#' This function is a wrapper around ProcessUniquePeptides with sensible defaults.
#'
#' @param df Dataframe
#' @param apply_filter If TRUE, filter non-unique peptides; if FALSE, only add the Unique column
#'
#' @return Modified dataframe
#'
#' @export
ApplyUniquePeptidesFilter <- function(df, apply_filter = TRUE) {
  
  # Create debug file in the current working directory
  debug_path <- file.path(getwd(), "data_with_unique_col.tsv")
  
  result <- ProcessUniquePeptides(
    df = df,
    save_debug = TRUE,
    debug_path = debug_path,
    filter_unique = apply_filter,
    peptide_col = "Precursor.Id",
    protein_col = "Protein.Group"
  )
  
  return(result)
}

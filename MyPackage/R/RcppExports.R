# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ReadBlock <- function(asciifname, start_row, numcols, numrows_in_block) {
    .Call('_Eagle_ReadBlock', asciifname, start_row, numcols, numrows_in_block)
}

ReadBlockBin <- function(binfname, start_row, numcols, numrows_in_block) {
    .Call('_Eagle_ReadBlockBin', binfname, start_row, numcols, numrows_in_block)
}

ReshapeM_rcpp <- function(fnameM, fnameMt, indxNA, dims) {
    .Call('_Eagle_ReshapeM_rcpp', fnameM, fnameMt, indxNA, dims)
}

calculateMMt_rcpp <- function(f_name_ascii, max_memory_in_Gbytes, num_cores, selected_loci, dims, quiet, message) {
    .Call('_Eagle_calculateMMt_rcpp', f_name_ascii, max_memory_in_Gbytes, num_cores, selected_loci, dims, quiet, message)
}

calculate_a_and_vara_batch_rcpp <- function(numreps, f_name_ascii, selected_loci, inv_MMt_sqrt, dim_reduced_vara, max_memory_in_Gbytes, dims, a, quiet, message) {
    .Call('_Eagle_calculate_a_and_vara_batch_rcpp', numreps, f_name_ascii, selected_loci, inv_MMt_sqrt, dim_reduced_vara, max_memory_in_Gbytes, dims, a, quiet, message)
}

calculate_a_and_vara_rcpp <- function(f_name_ascii, selected_loci, inv_MMt_sqrt, dim_reduced_vara, max_memory_in_Gbytes, dims, a, quiet, message) {
    .Call('_Eagle_calculate_a_and_vara_rcpp', f_name_ascii, selected_loci, inv_MMt_sqrt, dim_reduced_vara, max_memory_in_Gbytes, dims, a, quiet, message)
}

calculate_reduced_a_rcpp <- function(f_name_ascii, varG, P, y, max_memory_in_Gbytes, dims, selected_loci, quiet, message) {
    .Call('_Eagle_calculate_reduced_a_rcpp', f_name_ascii, varG, P, y, max_memory_in_Gbytes, dims, selected_loci, quiet, message)
}

createM_ASCII_rcpp <- function(f_name, f_name_ascii, type, AA, AB, BB, max_memory_in_Gbytes, dims, quiet, message, missing) {
    .Call('_Eagle_createM_ASCII_rcpp', f_name, f_name_ascii, type, AA, AB, BB, max_memory_in_Gbytes, dims, quiet, message, missing)
}

createM_BIN_rcpp <- function(f_name, f_name_bin, type, AA, AB, BB, max_memory_in_Gbytes, dims, quiet, message, missing) {
    .Call('_Eagle_createM_BIN_rcpp', f_name, f_name_bin, type, AA, AB, BB, max_memory_in_Gbytes, dims, quiet, message, missing)
}

createMt_ASCII_rcpp <- function(f_name, f_name_ascii, type, max_memory_in_Gbytes, dims, quiet, message) {
    invisible(.Call('_Eagle_createMt_ASCII_rcpp', f_name, f_name_ascii, type, max_memory_in_Gbytes, dims, quiet, message))
}

createMt_BIN_rcpp <- function(f_name_in, f_name_out, type, max_memory_in_Gbytes, dims, quiet, message) {
    invisible(.Call('_Eagle_createMt_BIN_rcpp', f_name_in, f_name_out, type, max_memory_in_Gbytes, dims, quiet, message))
}

extract_geno_rcpp <- function(f_name_ascii, max_memory_in_Gbytes, selected_locus, dims) {
    .Call('_Eagle_extract_geno_rcpp', f_name_ascii, max_memory_in_Gbytes, selected_locus, dims)
}

fasttimer <- function() {
    .Call('_Eagle_fasttimer')
}

getRowColumn <- function(fname) {
    .Call('_Eagle_getRowColumn', fname)
}


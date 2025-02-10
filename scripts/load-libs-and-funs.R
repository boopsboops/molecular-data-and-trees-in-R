#!/usr/bin/env Rscript

# cran
library("here")
library("glue")
library("tidyverse")
library("ape")
library("phangorn")
library("rentrez")
library("traits")
library("purrr")
library("ips")
library("castor")

# bioconductor
library("treeio")
library("tidytree")
library("ggtree")
library("DECIPHER")
library("msa")

# funs

# COLLAPSES HAPLOTYPES (FROM A DATAFRAME FORMAT TO A DATAFRAME FORMAT)
# TIDYVERSE VERSION
collapse_haplotypes <- function(df,lencol,seqcol) {
    adf <- df |> dplyr::arrange(dplyr::desc({{lencol}}))
    aseqs <- adf |> dplyr::pull({{seqcol}})
    reps <- aseqs |> purrr::map(\(x) which(stringr::str_detect(string=aseqs,pattern=x))[1]) |> unlist()
    ind <- unique(reps)
    adf.haps <- adf |> dplyr::slice(ind)
    dat.haps.ann <- adf.haps |> dplyr::mutate(nHaps=as.numeric(table(reps)))
    return(dat.haps.ann)
}

# BASE R VERSION
# COLLAPSES HAPLOTYPES (FROM A DATAFRAME FORMAT TO A DATAFRAME FORMAT)
hap_collapse_df <- function(df,lengthcol,nuccol){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    reps <- mapply(FUN=function(x) which(stringr::str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE)
    ind <- unique(reps)
    dat <- odf[ind,]
    dat[["nHaps"]] <- as.numeric(table(reps))
    return(dat)
}


# CONVERT TABULAR TO FASTA
tibble2dnabin <- function(df,seqcol,namecol) {
    seqs.list <- df |> dplyr::pull({{seqcol}}) |> stringr::str_split("")
    names(seqs.list) <- df |> dplyr::pull({{namecol}})
    seqs.list.dnabin <- seqs.list |> ape::as.DNAbin()
    return(seqs.list.dnabin)
}


# WRITE FASTA FROM A TIBBLE
write_tibble_fasta <- function(df,seqcol,namecol,file) {
    df2 <- df |> dplyr::mutate(fasta=glue::glue(">{.data[[namecol]]}\n{.data[[seqcol]]}"))
    file.create(file,overwrite=TRUE)
    df2 |> dplyr::pull(fasta) |> write(file=file,append=TRUE)
    writeLines(glue::glue("Fasta file written to '{file}'."))
}


# FUNCTION TO CONVERT DNABIN OBJECT TO TABULAR
dnabin2tibble <- function(dnas,namecol,seqcol) {
    dnas.list <- as.list(dnas)
    dnas.names <- names(dnas.list)
    dnas.char <- dnas.list |> purrr::map(\(x) paste(as.character(x),collapse=""))
    dnas.char.lower <- stringr::str_to_lower(dnas.char)
    dnas.clean <- stringr::str_replace_all(dnas.char.lower,"-|n|\\?","")
    dnas.tib <- tidyr::tibble({{namecol}}:=dnas.names,{{seqcol}}:=dnas.clean)
    return(dnas.tib)
}


# CLEAN GAPS OUT OF AN ALIGNMENT
clean_alignment_gaps <- function(matrix,maxgaps) {
    gaps.by.pos <- matrix |> as.character() |> apply(2,\(x) length(which(x=="-")))
    gaps.by.prop <- gaps.by.pos/dim(matrix)[1]
    gaps.keep <- which(gaps.by.prop<=maxgaps)
    matrix.reduced <- matrix[,gaps.keep]
    return(matrix.reduced)
}

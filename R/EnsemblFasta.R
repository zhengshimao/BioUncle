#' Read fasta file from Ensembl and get its information.
#'
#' Read fasta file from Ensembl and get its information.
#'
#' @title read_EnsemblFasta
#' @param fasta A fasta file from Ensembl, such as cDNA, pep, ncRNA.
#'
#' @return A data frame.
#' @export
read_EnsemblFasta <- function(fasta){
  # 1.  check the file
  two_line <- vroom::vroom_lines(file = fasta, n_max = 2,
                                 skip_empty_rows = TRUE)
  if( stringr::str_starts(two_line[1],pattern = ">") & !stringr::str_starts(two_line[2],pattern = ">") ){

    seq_char <- two_line[1] %>% stringr::str_remove(pattern = ">.*? ") %>% stringr::str_remove(pattern = " .*")

    check1 <- seq_char %in% c("pep","cdna","ncrna")
    check2 <- stringr::str_detect(two_line[1], pattern = "gene:|transcript:|gene_biotype:")

    if(check1 & check2 ){
      message("OK: the input file")
    }else{
      stop("The input file must be a fasta file of ncRNA/cDNA/peptide sequence from Ensembl!")
    }

  }else{
    stop("The input file must be a fasta file of ncRNA/cDNA/peptide sequence from Ensembl!")
  }


  # 2. read fasta
  fasta <- read_bio_fasta(fasta = fasta)


  # 3. extract info of fasta from Ensembl

  if(seq_char == "pep"){
    message("Note: Perhaps the input file is peptide sequence!")

    translation_id <- fasta$id %>% stringr::str_remove(pattern = " .*")

    transcript_id <- fasta$id %>%
      stringr::str_extract(pattern = "transcript:.*? ") %>%
      stringr::str_remove(pattern = "transcript:") %>%
      stringr::str_remove(pattern = " ")

  }else { # ncrna or cdna
    message("Note: Perhaps the input file is cDNA sequence!")
    transcript_id <- fasta$id %>% stringr::str_remove(pattern = " .*")
  }

  if(stringr::str_detect(two_line[1], pattern = "gene_symbol:")) # 对于模式生物的gene symbol
    gene_symbol <- fasta$id %>%
    stringr::str_extract(pattern = "gene_symbol:.*? ") %>%
    stringr::str_remove(pattern = "gene_symbol:") %>%
    stringr::str_remove(pattern = " ")

  position <- fasta$id %>% stringr::str_remove(pattern = ".*? .*? ") %>%
    stringr::str_remove(pattern = " .*") %>%
    stringr::str_remove(pattern = ".*?:.*?:")

  chr <- position %>% stringr::str_extract(pattern = ".*?:") %>%
    stringr::str_remove(pattern = ":")
  start <- position %>% stringr::str_extract(pattern = ".*?:.*?:") %>%
    stringr::str_remove(pattern = ".*?:") %>%
    stringr::str_remove(pattern = ":$")
  end <- position %>% stringr::str_remove(pattern = ".*?:.*?:") %>%
    stringr::str_remove(pattern = ":.*")
  strand <- position %>% stringr::str_remove(pattern = ".*:")

  gene_id <- fasta$id %>%
    stringr::str_extract(pattern = "gene:.*? ") %>%
    stringr::str_remove(pattern = "gene:") %>%
    stringr::str_remove(pattern = " ")

  transcript_biotype <- fasta$id %>%
    stringr::str_extract(pattern = "transcript_biotype:.*? ") %>%
    stringr::str_remove(pattern = "transcript_biotype:") %>%
    stringr::str_remove(pattern = " ")

  description <- fasta$id %>%
    stringr::str_extract(pattern = "description:.*") %>%
    stringr::str_remove(pattern = "description:")

  # 4. merge effective infomation
  if(seq_char == "pep"){

    if(stringr::str_detect(two_line[1], pattern = "gene_symbol:")){
      df_fasta <- data.frame(translation_id, gene_id,gene_symbol,transcript_id,chr, start, end, strand, transcript_biotype, description, sequence = fasta$sequence) %>%
        dplyr::mutate(len = stringr::str_length(sequence))
    }else{
      df_fasta <- data.frame(translation_id, gene_id,transcript_id,chr, start, end, strand, transcript_biotype, description, sequence = fasta$sequence) %>%
        dplyr::mutate(len = stringr::str_length(sequence))
    }

  }else{

    if(stringr::str_detect(two_line[1], pattern = "gene_symbol:")){
      df_fasta <- data.frame(transcript_id, chr, start, end, strand, gene_id,gene_symbol,transcript_biotype, description, sequence = fasta$sequence) %>%
        dplyr::mutate(len = stringr::str_length(sequence))
    }else{
      df_fasta <- data.frame(transcript_id, chr, start, end, strand, gene_id,transcript_biotype, description, sequence = fasta$sequence) %>%
        dplyr::mutate(len = stringr::str_length(sequence))
    }

  }

  return(df_fasta)
}

#' Get longest fasta from Ensembl cDNA/pep sequence.
#'
#' Get longest fasta and information from Ensembl cDNA/pep sequence file.
#'
#' @title get_LongestEnsemblFasta
#' @param fasta A fasta file of cDNA/pep sequence from Ensembl. Files ending in .gz, .bz2, .xz, or .zip will be automatically uncompressed.
#' @param source A reminder for you to use the Ensembl sequence.
#' @param longest_fa Write the longest fasta file. Files ending in .gz, .bz2, .xz, or .zip will be automatically compressed.
#' if `NULL`, it would be ignored.
#' @param seq_id Type of fasta IDs (default `gene_id`). `gene_id`, `translation_id`, `transcript_id`.
#' @param width An integer, (default `60`) number of bases of each line in the file.
#' @param longest_info Write the information of the longest fasta file. Files ending in .gz, .bz2, .xz, or .zip will be automatically compressed.
#' if `NULL`, it would be ignored.
#'
#' @return A list with two elements containing longest fasta data frame and its all information, `fa_longest` and `info_longest`.
#' @export
get_LongestEnsemblFasta <- function(fasta,
                                    source = "Ensembl",
                                    longest_fa = NULL,
                                    seq_id = "gene_id",
                                    width = 60,
                                    longest_info = NULL){
  # read fasta file from Ensembl cDNA or pep
  df_fasta <- read_EnsemblFasta(fasta = fasta)
  # get longest pep or cDNA
  message("Exacting longest fasta...")
  info_longest <- df_fasta %>% dplyr::group_by(gene_id) %>%
    dplyr::arrange(dplyr::desc(len), .by_group = T) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  # Confirm the type of sequence and check that the seq_id meets the specification
  seq_check <- df_fasta[1,"sequence"] %>% stringr::str_to_upper() %>% stringr::str_remove_all(pattern = "A|T|C|G|N")
  if(seq_check != ""){ # df_fasta: pep
    seq_id <- match.arg(seq_id, c("gene_id", "translation_id", "transcript_id"))
  }else{ # df_fasta: cdna or ncrna
    seq_id <- match.arg(seq_id, c("gene_id", "transcript_id"))
  }

  fa_longest <- info_longest %>% dplyr::select(all_of(seq_id), sequence)  # all_of: select the variables in the character vector

  if(!is.null(longest_fa))
    write_bio_fasta(fa_longest, file = longest_fa, width = width, append = F)

  if(!is.null(longest_info))
    vroom::vroom_write(info_longest, file = longest_info,
                       delim = "\t", eol = "\n", quote = "none")

  return(list(fa_longest = fa_longest,info_longest = info_longest))
}

#' Write the fasta file of longest sequence and its information.
#'
#' @title write_LongestEnsemblFasta
#' @param list_longest The result list of get_LongestEnsemblFasta function.
#' @param longest_fa Write the longest fasta file. Files ending in .gz, .bz2, .xz, or .zip will be automatically compressed.
#' if `NULL`, it would be ignored.
#' @param width An integer, (default `60`) number of bases of each line in the file.
#' @param longest_info Write the information of the longest fasta file. Files ending in .gz, .bz2, .xz, or .zip will be automatically compressed.
#' if `NULL`, it would be ignored.
#'
#' @export
write_LongestEnsemblFasta <- function(list_longest,
                                      longest_fa,
                                      width = 60,
                                      longest_info){
  write_bio_fasta(longest_list[["fa_longest"]], file = longest_fa, width = width, append = F)
  vroom::vroom_write(longest_list[["info_longest"]], file = longest_info,
                     delim = "\t", eol = "\n", quote = "none")
}

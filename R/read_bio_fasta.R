#' Create a data frame from a fasta file.
#'
#' Reads a file in fasta format and creates a data frame from it.
#'
#' @param fasta A fasta file.
#' @param simplified_ids A logical value(default:`FALSE`). `TRUE` for simplified ids of fasta sequence with characters between ‘>’ and first space, and `FALSE` for not simplified its ids.
#'
#' @return A data frame.
#' @export
read_bio_fasta <- function(fasta,
                           simplified_ids = FALSE
){
  fa <- vroom::vroom_lines(file = fasta, skip_empty_rows = TRUE, num_threads = 1)
  # a vector of fasta id
  id <- fa[stringr::str_detect(fa, "^>")] %>% stringr::str_remove(pattern = "^>")

  if(simplified_ids) id <- id %>% stringr::str_remove(pattern = " .*") # If the condition is TRUE, simplify ids of fasta sequence.

  # a vector of fasta sequence
  seq <- stringr::str_replace_all(fa,pattern = "^>.*",replacement = ">") %>% paste0(collapse = "") %>% stringr::str_split(pattern = ">")
  seq <- seq[[1]][-1]

  return(data.frame(id = id, sequence = seq))
}

#' Create a data frame from a fasta file.
#'
#' Reads a file in fasta format and creates a data frame from it.
#'
#' @param fasta A fasta file.
#'
#' @return A data frame.
#' @export
read_bio_fasta2 <- function(fasta
){
  fa <- vroom::vroom_lines(file = fasta, skip_empty_rows = TRUE)

  df_fa <- data.frame(name = character(0), seq = character(0))

  n=1 # Numbers of fasta sequences;row numbers of data frame.

  for (i in 1:length(fa)) {

    if(stringr::str_starts(fa[i],pattern = ">")){

      df_fa[n,1] <- fa[i] %>% stringr::str_remove(pattern = ">")
      n = n + 1

    }else{  # Not NA

      if(is.na(df_fa[n-1,2])){
        df_fa[n-1,2] <- fa[i]
      }else{
        df_fa[n-1,2] <- paste0(df_fa[n-1,2],fa[i])
      }

    }
  }

  return(df_fa)

}


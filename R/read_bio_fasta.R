#' Create a data frame from a fasta file.
#'
#' Reads a file in fasta format and creates a data frame from it.
#'
#' @param fasta A fasta file.
#'
#' @return A data frame.
#' @export
read_bio_fasta <- function(fasta
                           ){
  fa <- vroom::vroom_lines(file = fasta)

  df_fa <- data.frame(name = character(0), seq = character(0))

  n=1 # Numbers of fasta sequences;row numbers of data frame.

  for (i in 1:length(fa)) {

    if(stringr::str_starts(fa[i],pattern = ">")){

      df_fa[n,1] <- fa[i] %>% stringr::str_remove(pattern = ">")
      n = n + 1

    }else if(fa[i] != ""){  # Not NA

      if(is.na(df_fa[n-1,2])){
        df_fa[n-1,2] <- fa[i]
      }else{
        df_fa[n-1,2] <- paste0(df_fa[n-1,2],fa[i])
      }

    }
  }

  return(df_fa)

}
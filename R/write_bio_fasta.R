#' Write a data frame to a file in fasta format.
#'
#' `write_bio_fasta()`: write a data frame to a file in fasta format.
#'
#' @param data A data frame with 2 or 3 columns including `id`, `sequence` and `description`,separately.
#' @param file A path to a file of fasta sequence.
#' @param width An integer, (default `60`) numbers of bases of each line in the file.
#' @param append see \link[readr]{write_lines}. If `FALSE`, will overwrite existing file. If `TRUE`, will append to existing file. In both cases, if the file does not exist a new file is created.
#' @param num_threads see \link[readr]{write_lines}. The number of processing threads to use for initial parsing and lazy reading of data. If your data contains newlines within fields the parser should automatically detect this and fall back to using one thread only. However if you know your file has newlines within quoted fields it is safest to set `num_threads = 1` explicitly.
#'
#' @return A file in fasta format.
#' @export
#'
#' @examples
#' \dontrun{
#' data(pst_hrp)
#' write_bio_fasta(data = pst_hrp, file = "pst_hrp.fasta")
#' }
write_bio_fasta <- function(data,
                            file,
                            width = 60,
                            append = FALSE,
                            num_threads = readr::readr_threads()
                            ){

  # check data
  stopifnot(is.data.frame(data))
  n_col <- ncol(data)
  if(n_col != 2 & n_col != 3){
    stop("`data`: the input must be a data frame with 2 or 3 columns!")
  }
  witdh = as.integer(width)
  # prepare the list for `write_lines()`
  fa <- list()
  for (i in seq(nrow(data))){

    n = length(fa) + 1 #向fa list 添加到第n个元素，也就是fasta序列文件的行数了。

    if(n_col == 2){
      id <- data[i,1] %>% stringr::str_replace(pattern = "^", ">") %>% as.character()
      fa[[n]] <- id
    } else {
      id <- data[i,1] %>% stringr::str_replace(pattern = "^", ">") %>% as.character()
      description <- data[i,3] %>% as.character()
      fa[[n]] <- paste(id,description,sep = " ")
    }

    seq <- data[i,2] %>% as.character()

    for (x in seq(ceiling(stringr::str_length(seq)/width))) {
      fa[[n+1]] <- stringr::str_sub(seq,start = ((x-1)*width+1), end = x*width)
      n=length(fa)
    }
  }

  # write out `fa` list
  readr::write_lines(fa,file = file, sep = "\n", append = FALSE, num_threads = readr::readr_threads()) #要求为list或向量。若为list，写出的list内容，不含name
  invisible(gc())
}

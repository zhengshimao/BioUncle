#' Convert to single cell according to the same content.
#'
#' A function opposite to \link[tidyr]{separate_rows}.
#'
#' @param df A data frame with two columns.
#' @param sep A separator character.
#' @param rename A vector of new column names of length 2.
#'
#' @return A data frame with two character vectors.
#' @export
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' df <- iris %>% select(5,3)
#' str(df)
#'
#' new_df <- get_one2mul(df)
#' str(new_df)
#' }
get_one2mul <- function(df,
                        sep = "/",
                        rename = NULL
                        ){
  stopifnot(is.data.frame(df)) # df
  stopifnot(ncol(df) == 2) # df: 2 columns
  if(!is.null(rename)){stopifnot(length(rename) == 2)}

  df_name <- names(df)
  names(df) <- c("A","B")
  df$A <- df$A %>% as.character()
  df$B <- df$B %>% as.character()
  # a new data frame with two columns
  res <- data.frame(A = character(0), B = character(0)) #two character vecttors.

  for (i in 1:nrow(df)) {

    if(df[i,1] %in% res[,1]){
      position_nrow <- which(res$A %in% df[i,1])
      res[position_nrow,2] <- paste(res[position_nrow,2],df[i,2], sep = "/")
    }else{
      # a new content
      res[nrow(res)+1,1] <- df[i,1]
      res[nrow(res),2] <- df[i,2] #don't add the number of row again.
    }

  }
  if(is.null(rename)){
    names(res) <- df_name
  }else{
    names(res) <- rename
  }

  return(res)

}

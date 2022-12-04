#' Convert to single cell according to the same content.
#'
#' A function opposite to \link[tidyr]{separate_rows}.
#'
#' @param df A data frame with two columns.
#' @param sep A separator character.
#'
#' @return A data frame with two character vectors.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- iris %>% select(5,3)
#' str(df)
#' new_df <- get_one2mul(df)
#'
#' str(new_df)
#' }
get_one2mul <- function(df,
                        sep = "/"
                        ){
  stopifnot(is.data.frame(df)) # df

  stopifnot(ncol(df) == 2) #能力暂时只能做2列数据的，凑合用吧

  df_name <- names(df)
  names(df) <- c("A","B")

  df$A <- df$A %>% as.character()
  res <- data.frame(A = character(0), B = character(0)) #结果每列会变为字符串向量。可以根据输入df判断写活。

  for (i in 1:nrow(df)) {

    if(df[i,1] %in% res[,1]){
      position_nrow <- which(res$A %in% df[i,1])
      res[position_nrow,2] <- paste(res[position_nrow,2],df[i,2], sep = "/")
    }else{
      # gene 首次出现
      res[nrow(res)+1,1] <- df[i,1]
      res[nrow(res),2] <- df[i,2] #可不要再给行数＋1了
    }

  }

  names(res) <- df_name
  return(res)

}

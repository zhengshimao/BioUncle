#' Convert FPKM to TPM.
#'
#' Convert FPKM matrix from RNA-seq to TPM matrix.
#'
#' @param FPKM A FPKM matrix or data frame of RNA-seq expression.
#'
#' @return A data frame.
#' @export
fpkm2tpm <- function(FPKM){

  if(!(is.matrix(fpkm) | is.data.frame(fpkm)) & all(fpkm>=0))
    cli::cli_abort(c("x" = "FPKM must be a matrix or data frame!","x" = "Each value in FPKM must >=0."))

  # 转换
  FPKM2TPM <- function(fpkm){
      exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }

  TPM <- apply(FPKM,2,FPKM2TPM) %>% as.data.frame()
  # tryCatch({TPM <- apply(FPKM,2,FPKM2TPM) %>% as.data.frame()},
  #          error = function(e){
  #            cli::cli_abort(c("i" = "NOTE: FPKM matrix could contain character columns, but not numeric columns!"))
  #          },
  #          finally = {
  #            TPM
  #          }
  #          )

  cli::cli_alert_success("There {?is/are} {ncol(TPM)} sample{?s} in total!")
  return(TPM)
}







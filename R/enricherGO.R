#' A more user-friendly enrichment analyzer based on enricher function.
#'
#' The brother of enricher and enrichGO function.
#'
#' @title enricherGO
#' @param gene See \link[clusterProfiler]{enricher}.
#' @param go2ont a data frame with "go_id’ and "Ontology". If `NULL`, it will be downloaded every time. It is recommended to save the `go2ont` with the `get_GO2Ontology_table()` function first.
#' @param ont GO Ontology: "BP", "MF", "CC", "ALL".
#' @param pvalueCutoff See \link[clusterProfiler]{enricher}.
#' @param qvalueCutoff See \link[clusterProfiler]{enricher}.
#' @param pAdjustMethod See \link[clusterProfiler]{enricher}.
#' @param universe See \link[clusterProfiler]{enricher}.
#' @param minGSSize See \link[clusterProfiler]{enricher}.
#' @param maxGSSize See \link[clusterProfiler]{enricher}.
#' @param TERM2GENE A data frame with GO IDs and gene IDs. See \link[clusterProfiler]{enricher}.
#' @param TERM2NAME A data frame with GO IDs and GO descriptions. See \link[clusterProfiler]{enricher}.
#'
#' @return A object of enrichResult
#' @export
enricherGO <- function(gene,
                       go2ont = NULL,
                       ont = "ALL",
                       pvalueCutoff,
                       qvalueCutoff = 0.2,
                       pAdjustMethod = "BH",
                       universe,
                       minGSSize = 10,
                       maxGSSize = 500,
                       TERM2GENE,
                       TERM2NAME = NA
                       ){
  # go2ont
  if(is.null(go2ont)){
    message("Downloading……\nGO2Ontology table.")
    go2ont = clusterProfiler:::get_GO2Ontology_table() # devtools::check()  NOTE!
  }else{
    go2ont =  as.data.frame(go2ont)
    names(go2ont) <- c("go_id","Ontology")

    if(length(unique(go2ont[,"Ontology"])) == 1){
      ont <- unique(go2ont[,"Ontology"])
    }

    message("GO2Ontology: Ready!")
  }
  ont <- stringr::str_to_upper(ont)
  ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

  if(ont != "ALL")
    go2ont <- go2ont %>% dplyr::filter(Ontology == ont)

  # enricher
  go_result <-  clusterProfiler::enricher(
    gene = gene,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME)

  if(ont == "ALL"){
    go_result@ontology <- "GOALL"
  }else{
    go_result@ontology <- ont
  }

  # change the order of colnames
  go_result@result <- go_result@result %>% dplyr::left_join(go2ont,by = c("ID" = "go_id")) %>% dplyr::select(10,1,2,9,3:8) %>% dplyr::arrange(pvalue) %>% tidyr::drop_na()
  # rename,mkae it same to the result of enrichGO function.
  colnames(go_result@result)[1] <- "ONTOLOGY"

  return(go_result)
}


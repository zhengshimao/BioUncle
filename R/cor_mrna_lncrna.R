#' Calculate the correlation between the two expression matrix.
#'
#' cor_mrna_lncrna(): Calculate the correlation between the two expression matrix.
#'
#' @param filter_log2_mrna Expression matrix of mRNA transformed by log2.
#' @param filter_log2_lncrna Expression matrix of lncRNA transformed by log2.
#' @param cor_method The method of correlation. See \link[stats]{cor}.
#' @param padj_method Correction method. See \link[stats]{p.adjust}
#' @param digits See \link[base]{round}.
#' @param check_mat Logical. Check the expression matrix or not.
#' @param pstar Logical. paste stars
#'
#' @return A list
#' @export
#'
cor_mrna_lncrna <- function(filter_log2_mrna,
                            filter_log2_lncrna,
                            cor_method = "pearson",
                            padj_method = "BH",
                            digits = 2,
                            check_mat = TRUE,
                            pstar = FALSE){

  # 释放内存
  invisible(gc())

  # 检查参数
  # 就检查两个矩阵吧
  digits <- as.integer(digits)
  stopifnot(is.integer(as.integer(digits)) & as.integer(digits) > 0)

  if(check_mat){
    if(!all(filter_log2_mrna >= 0)){
      stop("The matrix of filter_log2_mrna must be nonnegative number!\n",
           "There are ",sum(filter_log2_mrna < 0), " numbers less than 0!")
    }

    if(!all(filter_log2_lncrna >= 0)){
      stop("The matrix of filter_log2_lncrna must be nonnegative number!\n",
           "There are ",sum(filter_log2_lncrna < 0), " numbers less than 0!")
    }

    if(any(filter_log2_mrna >= 20)){
      warning("The matrix of filter_log2_mrna may transformed using 'log2(x+1)'!\n",
              "There are ",sum(filter_log2_mrna >= 20), " numbers greater than 20!")

    }

    if(any(filter_log2_lncrna >= 20)){
      warning("The matrix of filter_log2_lncrna may be transformed using 'log2(x+1)'!\n",
              "There are ",sum(filter_log2_lncrna >= 20), " numbers greater than 20!")
    }
  }

  # 主体
  # suppressPackageStartupMessages(library(magrittr))
  # 检查与调整两个输入矩阵的样本名顺序，并调整为一致
  if(!identical(colnames(filter_log2_lncrna),colnames(filter_log2_mrna))){
    filter_log2_mrna <- filter_log2_mrna[,match(colnames(filter_log2_lncrna),colnames(filter_log2_mrna))]
  }else{
    message("Samplenames: consistent!")
  }
  # 相关性矩阵计算
  cor_r <- WGCNA::cor(t(filter_log2_mrna),t(filter_log2_lncrna), method = cor_method)
  # 未矫正p值计算
  wgcna_p <- WGCNA::corPvalueStudent(cor = cor_r,nSamples = ncol(filter_log2_mrna))
  # 矫正p值计算
  padj <- stats::p.adjust(wgcna_p, method = padj_method) %>% base::matrix(nrow = nrow(wgcna_p))
  rownames(padj) <- rownames(cor_r)
  colnames(padj) <- colnames(cor_r)
  # 生成部分结果res，设为list结构
  res <- list()
  res$r <- cor_r
  res$p <- wgcna_p
  res$padj <- padj
  res$padj_method <- padj_method

  # 释放内存
  invisible(gc())

  # 计算r+stars，自定义函数
  pstars <- function(r, p, digits = digits){
    r <- round(r,digits = digits) # matrix

    #维度检查
    if(!all(dim(r) == dim(p))){
      stop("Dimensions of r and p matrix are not identical!")
    }
    # 检查行列名
    if(!all(rownames(r) %in% rownames(p))){
      stop("Rownames of r and p matrix are not identical!")
    }

    if(!all(colnames(r) %in% colnames(p))){
      stop("Colnames of r and p matrix are not identical!")
    }

    # 检查r与p行列名排序。使得行列顺序一致
    p_row <- identical(rownames(r),rownames(p))
    p_col <- identical(colnames(r),colnames(p))
    if(!(p_row & p_col)){
      p <- p[match(rownames(r),rownames(p)),match(colnames(r),colnames(p))]
    }
    # 计算
    message("get stars!")
    stars0 <- ifelse(p < 0.001,"***",
                     ifelse(p < 0.01,"**",
                            ifelse(p < 0.05,"*","")
                     )
    ) # matrix

    stars <- paste0(r,stars0) %>% matrix(nrow = nrow(r))
    rownames(stars) <- rownames(r)
    colnames(stars) <- colnames(r)

    return(stars)
  }

  # 生成部分结果res，设为list结构
  if(pstar){
    res$p_stars <- pstars(res$r,res$p, digits = digits)
    res$padj_stars <- pstars(res$r,res$padj, digits = digits)
  }else{
    res$p_stars <- NULL
    res$padj_stars <- NULL
  }

  # 释放内存
  invisible(gc())

  return(res)

}

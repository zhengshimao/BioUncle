#' Write the result of `cor_mrna_lncrna()` to disk.
#'
#' write_cor_p(): write the result of `cor_mrna_lncrna()` to disk.
#'
#' @param cor_p A list from \link[BioUncle]{cor_mrna_lncrna}.
#' @param r_cutoff r cutoff.
#' @param padj_cutoff Padj cutoff.
#' @param output_path The path of results files.
#' @param prefix The prefix of results files.
#' @param num_threads Number of threads, see\link[vroom]{vroom_write}.
#' @param gzip Logical, compress the output or not.
#' @param get_WiderData Logical, write wide data frame or not.
#'
#' @return Result files.
#' @export
write_cor_p <- function(cor_p,
                        r_cutoff = 0.9,
                        padj_cutoff = 0.05,
                        output_path = "./",
                        prefix = "my",
                        num_threads = 5,
                        gzip = TRUE,
                        get_WiderData = FALSE
                        ){
  #####主体
  # suppressPackageStartupMessages(library(data.table))
  # suppressPackageStartupMessages(library(tibble)) #行名转为列
  # suppressPackageStartupMessages(library(R.utils)) #压缩文件
  # suppressPackageStartupMessages(library(tidyr)) #宽数据变长数据
  # suppressPackageStartupMessages(library(dplyr))
  # suppressPackageStartupMessages(library(vroom))

  # 释放内存
  invisible(gc())

  # chek input
  if(!dir.exists(output_path)){
    message("DONE: create ", output_path,"\n")
    dir.create(output_path)
  }else{
    message("OK: ",output_path," is ready!\n")
  }

  ###### 处理为长数据集
  message("Converting r matrix !\n")
  long_r <- as.data.frame(cor_p$r) %>% tibble::rownames_to_column(var = "mRNA") %>% tidyr::pivot_longer(2:ncol(cor_p$r),names_to = "lncRNA",values_to = "r")
  message("Converting p matrix !\n")
  long_p <- as.data.frame(cor_p$p) %>% tibble::rownames_to_column(var = "mRNA") %>% tidyr::pivot_longer(2:ncol(cor_p$p),names_to = "lncRNA",values_to = "p")
  message("Converting padj matrix !\n")
  long_padj <- as.data.frame(cor_p$padj) %>% tibble::rownames_to_column(var = "mRNA") %>% tidyr::pivot_longer(2:ncol(cor_p$padj),names_to = "lncRNA",values_to = "padj")
  message("Combining r+p+padj matrix !\n")
  long_r_p <- long_r %>% dplyr::full_join(long_p, by=c("mRNA","lncRNA")) %>% dplyr::full_join(long_padj,by=c("mRNA","lncRNA")) %>% dplyr::arrange(dplyr::desc(abs(r)),p)
  message("Filtering r+p+padj matrix using r_cutoff = ",r_cutoff," and padj_cutoff = ",padj_cutoff,"!\n")
  long_rcut_padjcut <- long_r_p %>% dplyr::filter(r >= r_cutoff, padj <= padj_cutoff)
  invisible(gc()) #清理内存。有执行效果，还能压制信息

  ###### 写出数据
  if(get_WiderData){  ###### 是否写出宽数据部分

    r_file <- paste0(output_path,"/",prefix,"_r.tsv")
    if(gzip) r_file <- paste0(r_file,".gz")
    message("Writing ",r_file," to disk!\n")
    cor_p$r <- cor_p$r %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_id")
    # data.table::fwrite(cor_p$r,file = r_file,sep = "\t", na = "", quote = "auto",row.names = F, col.names = T,nThread = nThread)
    vroom::vroom_write(cor_p$r,file = r_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

    # if(gzip){
    #   if(file.exists(paste0(r_file,".gz"))) file.remove(paste0(r_file,".gz"))
    #   R.utils::gzip(r_file,remove = T)
    # }


    p_file <- paste0(output_path,"/",prefix,"_p.tsv")
    if(gzip) p_file <- paste0(p_file,".gz")
    message("Writing ",p_file," to disk!\n")
    cor_p$p <- cor_p$p %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_id")
    # data.table::fwrite(as.data.frame(cor_p$p),file = p_file,sep = "\t", na = "", quote = "auto", row.names = F, col.names = T,nThread = nThread)
    vroom::vroom_write(cor_p$p,file = p_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

    # if(gzip){
    #   if(file.exists(paste0(p_file,".gz"))) file.remove(paste0(p_file,".gz"))
    #   R.utils::gzip(p_file,remove = T)
    # }

    padj_file <- paste0(output_path,"/",prefix,"_padj.tsv")
    if(gzip) padj_file <- paste0(padj_file,".gz")
    message("Writing ",padj_file," to disk!\n")
    cor_p$padj <- cor_p$padj %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_id")
    # data.table::fwrite(as.data.frame(cor_p$padj),file = padj_file,sep = "\t", na = "", quote = "auto", row.names = F, col.names = T,nThread = nThread)
    vroom::vroom_write(cor_p$padj,file = padj_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

    # if(gzip){
    #   if(file.exists(paste0(padj_file,".gz"))) file.remove(paste0(padj_file,".gz"))
    #   R.utils::gzip(padj_file,remove = T)
    # }

    if(!is.null(cor_p$p_stars)){
      p_stars_file <- paste0(output_path,"/",prefix,"_p_stars.tsv")
      if(gzip) p_stars_file <- paste0(p_stars_file,".gz")
      message("Writing ",p_stars_file," to disk!\n")
      cor_p$p_stars <- cor_p$p_stars %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_id")
      # data.table::fwrite(as.data.frame(cor_p$p_stars),file = p_stars_file,sep = "\t", na = "", quote = "auto", row.names = F, col.names = T,nThread = nThread)
      vroom::vroom_write(cor_p$p_stars,file = p_stars_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

      # if(gzip){
      #   if(file.exists(paste0(p_stars_file,".gz"))) file.remove(paste0(p_stars_file,".gz"))
      #   R.utils::gzip(p_stars_file,remove = T)
      # }
    }

    if(!is.null(cor_p$padj_stars)){
      padj_stars_file <- paste0(output_path,"/",prefix,"_padj_stars.tsv")
      if(gzip) padj_stars_file <- paste0(padj_stars_file,".gz")
      message("Writing ",padj_stars_file," to disk!\n")
      cor_p$padj_stars <- cor_p$padj_stars %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_id")
      # data.table::fwrite(as.data.frame(cor_p$padj_stars),file = padj_stars_file,sep = "\t", na = "", quote = "auto",row.names = F, col.names = T,nThread = nThread)
      vroom::vroom_write(cor_p$padj_stars,file = padj_stars_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

      # if(gzip){
      #   if(file.exists(paste0(padj_stars_file,".gz"))) file.remove(paste0(padj_stars_file,".gz"))
      #   R.utils::gzip(padj_stars_file,remove = T)
      # }
    }
  }

  long_r_p_file <- paste0(output_path,"/",prefix,"_long_r_p.tsv")
  if(gzip) long_r_p_file <- paste0(long_r_p_file,".gz")
  message("Writing ",long_r_p_file," to disk!\n")
  # data.table::fwrite(long_r_p,file = long_r_p_file,sep = "\t", na = "", quote = "auto", row.names = F, col.names = T,nThread = nThread)
  vroom::vroom_write(long_r_p,file = long_r_p_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

  # if(gzip){
  #   if(file.exists(paste0(long_r_p_file,".gz"))) file.remove(paste0(long_r_p_file,".gz"))
  #   R.utils::gzip(long_r_p_file,remove = T)
  # }

  long_rcut_padjcut_file <- paste0(output_path,"/",prefix,"_long_r",r_cutoff,"_padj",padj_cutoff,".tsv")
  if(gzip) long_rcut_padjcut_file <- paste0(long_rcut_padjcut_file,".gz")
  message("Writing ",long_rcut_padjcut_file," to disk!\n")
  #data.table::fwrite(long_rcut_padjcut,file = long_rcut_padjcut_file,sep = "\t", na = "", quote = "auto", row.names = F, col.names = T,nThread = nThread)
  vroom::vroom_write(long_rcut_padjcut,file = long_rcut_padjcut_file, delim = "\t", na = "", quote = "none", num_threads = num_threads)

  # if(gzip){
  #   if(file.exists(paste0(long_rcut_padjcut_file,".gz"))) file.remove(paste0(long_rcut_padjcut_file,".gz"))
  #   R.utils::gzip(long_rcut_padjcut_file,remove = T)
  # }

  # 释放内存
  invisible(gc())
}

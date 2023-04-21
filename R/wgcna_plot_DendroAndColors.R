#' Plot dendrogram and mudule colors.
#'
#' Plot dendrogram and mudule colors using ggtree. This function is similar to WGCNA::plotDendroAndColors() and
#' It origines from `ggtreeDendro::plot_wgcna()`.
#'
#' @param x Output of `WGCNA::blockwiseModules()`.
#' @param abline_height Abline height. If it is `NULL`, abline would not be added.
#' @param abline_color The color of abline.
#' @param n_guide_lines Number of guide lines. If it is `NULL`, guide lines would not be added.
#' @param gudide_linetype Specifiy the type of line segment. Learn more about setting these aesthetics in vignette("ggplot2-specs").
#' @param line_interval The distance between tips and guide lines.
#' @param line_lowest_height Specifiy the lowest height of all guide lines. "auto" or a numberic value.
#' @param offset_tree_heatmap The distance between tree object and heatmap.
#' @param heamap_height Heatmap height. "auto" or a numberic value.
#' @param group_labels Label for the heatmap.
#' @param group_labels_position Position of `group_labels`. "left" or "right".
#'
#' @return A ggtree object.
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_segment
#' @importFrom ggtree gheatmap
#' @importFrom ggtree ggtree
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggh4x guide_axis_truncated
#' @export
wgcna_plot_DendroAndColors <- function(x,
                                       # abline
                                       abline_height = NULL,
                                       abline_color = "red",
                                       # guide lines
                                       n_guide_lines = 50,
                                       gudide_linetype = "dotted",
                                       line_interval = 0,
                                       line_lowest_height = "auto",
                                       # heatmap
                                       offset_tree_heatmap = 0,
                                       heamap_height = "auto",
                                       group_labels = "Mudule colors",
                                       group_labels_position = "left"
                                       ){
  # 1.plot basic dendrogram

  p <- ggtree(x$dendrograms[[1]], layout = "dendrogram", ladderize = FALSE, size = 0.2)+
    labs(x = "Height", title = "Cluster Dendrogram")+
    theme(axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text(),
          text = element_text(family = "serif", face = "bold"),
          plot.title = element_text(hjust = 0.5))+
    # scale_y_continuous(breaks = lables_y)+ # Why doesn't it work?
    guides(y = "axis_truncated")

  # 2.add abline

  if(is.numeric(abline_height) && length(abline_height) == 1)
    p <- p + geom_vline(xintercept = -abline_height,
                        color = abline_color) ## Note: It is "vline" not "hline".

  # 3.add giude lines

  if(is.numeric(n_guide_lines) && length(n_guide_lines) == 1){

    if(n_guide_lines != as.integer(n_guide_lines))
      cli::cli_warn(c("`n_guide_lines` is forced to {as.integer(n_guide_lines)}!",
                      "i" = "`n_guide_lines` must be an integer!"))

    ## extract the positons of tips from p object and sort it.
    gene_label <- p$data %>% dplyr::filter(isTip == TRUE) # extract coordinate of tips from p object.

    gene_label <- gene_label[x$dendrograms[[1]]$order,] # sort labels acording to x$dendrograms[[1]]$order.

    ## ascertain the positions of guide lines
    # lineSpacing = floor(length(x$dendrograms[[1]]$order)/(n_guide_lines + 1));
    lineSpacing = floor(length(x$dendrograms[[1]]$order)/(n_guide_lines - 1));
    positions = c(1, (1:(n_guide_lines-2))* lineSpacing, length(x$dendrograms[[1]]$order));

    if(line_lowest_height == "auto"){

      # line_lowest_height <- ((min(x$dendrograms[[1]]$height) - 0.05) %>% round(.,digits = 1)) - 0.1
      # line_lowest_height <- round(min(x$dendrograms[[1]]$height) - 0.05, digits = 1) - 0.1

      height_min <- min(x$dendrograms[[1]]$height)
      height_max <- max(x$dendrograms[[1]]$height)
      line_lowest_height <- height_min - (height_max - height_min)*0.2
    }
    else if(is.numeric(line_lowest_height) && length(line_lowest_height) == 1)
      line_lowest_height <- line_lowest_height
    else
      cli::cli_abort(message = c("x" = "`line_lowest_height` must be \"auto\" or a numeric value!"))

    guide_lines <- gene_label[positions,] %>% dplyr::select(x,y) %>%  #Coordinate data frame of tips.
      dplyr::mutate(guide_len = x - (-line_lowest_height) - line_interval)


    ## add guide lines
    p <- p + geom_segment(data = guide_lines,
                          mapping = aes(x = x, xend = x - guide_len, y = y, yend = y),
                          linetype  = gudide_linetype)
    # p # Warning message:  Removed 50 rows containing missing values (`geom_segment()`).
  }else if(!is.null(n_guide_lines) ){
    cli::cli_abort(c("x" = "`n_guide_lines` must be an integer!"))
  }

  # 4. add heatmap

  ## 4.1 get Module colors
  colors <- x$colors

  if(is.numeric(colors))
    colors <- WGCNA::labels2colors(colors)

  ## 4.2 make dataframe of heatmap
  heat_label <- c("module")
  names(heat_label) <- group_labels

  d <- data.frame(module=as.character(colors),
                  row.names=seq_along(colors)) %>%
    dplyr::rename( all_of(heat_label) )

  ## 4.3 define position of heatmap group_labels
  group_labels_position <- match.arg(group_labels_position, choices = c("left","right"))

  if(group_labels_position == "left"){
    colnames_position = "bottom";
    hjust = 1
  }else{
    colnames_position = "top";
    hjust = 0
  }

  if(heamap_height == "auto")
    heamap_height <- (height_max - height_min)*0.25
  else
    heamap_height <- heamap_height
  # 5. plot complete tree and heatmap
  gheatmap(p, d, color=NA,
           width = heamap_height,
           offset = offset_tree_heatmap,
           # colnames = TRUE,colnames_angle = 0,
           colnames_position = colnames_position,
           hjust = hjust) +
    scale_fill_identity() +
    labs(caption = as.character(as.expression(x$dendrograms[[1]]$call)))

}

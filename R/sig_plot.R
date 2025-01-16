#' Plotting function for mmsig package
#'
#' @param sig Object generated through \code{mm_fit_signatures}
#' @param plot.type Type of barplot. Options are 'stacked' (default) or '96'.
#' @param sig_order Order of signatures in \code{stacked} barplot, from top to bottom.
#' @param plot.names boolean whether or not to plot sample names on the x axis. Default = \code{FALSE}
#' @param col.width Widht for col bars. Default = 0.9
#'
#' @return relative contribution of mutational signatures in each sample
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 labs
#' @export
#'
plot_signatures = function(sig, plot.type = "stacked",
                           sig_order = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9", "SBS18", "SBS-MM1", "SBS35"),
                           plot.names = FALSE, col.width = 0.9){

  if(plot.type == "stacked"){
    sigPlot <- sig$estimate %>%
    rownames_to_column(var = "sample") %>%
    dplyr::select(-mutations) %>%
    melt(id.var = "sample", variable.name = "SBS", value.name = "prop") %>%
    mutate(SBS = factor(SBS, levels = sig_order)) %>%
    ggplot(aes(sample, prop, fill = SBS)) +
      geom_col(width = col.width)+
      scale_fill_sigs()+ scale_y_continuous(expand = c(0,0))+
      labs(x = "Sample", y = "Relative contribution", fill = "Signature")+
      theme_bw()+ theme(text = element_text(size = 10, color = "black"),
            axis.text = element_text(color = "black"),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            axis.title.x = element_blank())

      if(plot.names){
        sigPlot +
          theme(axis.text.x = rotatedAxisElementText(90, "top"))
      } else{
        sigPlot +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
      }
    }else if(plot.type == "96"){
      sig$mutmatrix %>% 
      mutate(mut = rownames(.), mutgroup = gsub(".*.\\[|\\].*.", "", rownames(.))) %>% 
      reshape2::melt() %>% 
      ggplot(aes(x = mut, y = value, fill = mutgroup)) + geom_col() +
        facet_grid(variable ~ mutgroup, scales = "free_x") + theme_bw() + 
        scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank(), legend.position = "none")
    }else{
      stop("ERROR: Please, indicate a correct plot type (ie, 'stacked' or '96' (as character, no number))")
    }
}

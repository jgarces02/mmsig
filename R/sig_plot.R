#' Plotting function for mmsig package
#'
#' @param sig Object generated through \code{mm_fit_signatures}
#' @param plot.type Type of barplot. Options are 'stacked' (default) or '96', or 'boot' for plotting boostrapped results.
#' @param plot.sample Indicate the same of the sample(s) to plot. Default = \code{"all"} (ie, plot all samples).
#' @param sig_order Order of signatures in \code{stacked} barplot, from top to bottom.
#' @param plot.names boolean whether or not to plot sample names on the x axis. Default = \code{FALSE}
#' @param col.width Widht for col bars. Default = 0.9
#'
#' @return relative contribution of mutational signatures in each sample
#' @importFrom dplyr mutate
#' @importFrom ggplot2
#' @export
#'
plot_signatures = function(sig, plot.type = "stacked", plot.sample = "all",
                           sig_order = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9", "SBS18", "SBS-MM1", "SBS35"),
                           plot.names = FALSE, col.width = 0.9){
  
  mysamples <- ifelse(plot.sample == "all", rownames(sig$estimate$sample), plot.sample)
  
  if(plot.type == "stacked"){
    sigPlot <- sig$estimate %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select(-mutations) %>% subset(sample %in% mysamples) %>%
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
      reshape2::melt() %>% subset(variable %in% mysamples) %>%
      ggplot(aes(x = mut, y = value, fill = mutgroup)) + geom_col() + geom_hline(yintercept = 0) +
        facet_grid(variable ~ mutgroup, scales = "free") + ggplot2::coord_cartesian(expand = F) +
        theme_bw() + labs(x = "SBS signature (tri-nt change)", y = "# mutations") +
        scale_fill_manual(values = c("lightskyblue", "black", "firebrick2","gray88","darkolivegreen3", "lightpink1")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank(), panel.border = element_blank(),
              legend.position = "none", axis.text.y = element_text(size = 12),
              panel.spacing = unit(0,'lines'), strip.background = element_blank(), strip.text.x = element_blank())
    
  }else if(plot.type == "boot"){
    ggplot(sig$bootstrap, aes(signature, estimate, fill = signature))+
      geom_bar(position="dodge", stat="identity")+
      geom_errorbar(data = sig$bootstrap, mapping = aes(x = signature, ymin = CI025, ymax = CI975))+
      facet_grid(~sample)+
      theme_bw() + theme(text = element_text(size = 12),
                         axis.text.x = rotatedAxisElementText(90, 'top'),
                         axis.title.x = element_blank(),
                         strip.background = element_blank(),
                         legend.position = 'none')+
      scale_fill_sigs()+ labs(y = 'Relative contribution') +
      ggplot2::coord_cartesian(ylim = c(0,1), expand = F)
  }else{
    stop("ERROR: Please, indicate a correct plot type (ie, 'stacked', '96' (as character, no number), or 'boot')")
  }
}

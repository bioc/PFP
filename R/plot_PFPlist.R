globalVariables(c("net", "sim", "refnet_index"))

# library(tidyr)
# library(plyr)

#' Plot multiple PFPs.
#'
#' Function for visualization multiple PFPs.
#'
#' @param object, \code{PFP} a list of PFP.
#' @param l_size, line size of plot, default is 0.5.
#' @aliases plot_PFPlist
#' @seealso \code{\link{PFP-class}}
#' @return plot the PFP list
#' @examples
#' data(PFP_test1)
#' pfp_list <- list(a=PFP_test1)
#' plot_PFPlist(pfp_list)
#' @export
plot_PFPlist <- function(object,
                         l_size = 0.5){
    if (is.list(object))
    {
        pfp_num <- length(object)
        pfp_score <- llply(object,
                           function(object)pathways_score(object)[["PFP_score"]]) %>%
            do.call(what = cbind) %>%
            as.data.frame
        network_num <- nrow(pfp_score)
        sim_df <- data.frame(pfp_score,refnet_index = seq(1):network_num) %>%
            gather(net,sim,-refnet_index)
        p <- ggplot(data = sim_df,
                    aes(x = refnet_index, y =sim)) +
            geom_line(aes(color = net),size = l_size) +
            ylab('Pathway fingerprint') + xlab('Index of basic Pathways') +
            guides(color = guide_legend("query network"))
        print(p)
    }
    else
        stop('object must be a PFP list while visualization of multiple pfps')
}

#' Title
#'
#' @param topoFile
#' @param matOut
#' @param plotOut
#' @param writeOut
#' @param nodeLabelSize
#'
#' @return
#' @export
#'
#' @examples
correlationMatBool <- function(topoFile,
                               matOut = F,
                               plotOut = F,
                               writeOut = T,
                               nodeLabelSize = NULL)
{
    print(topoFile)
    net <- topoFile %>% str_remove(".topo")
    if (file.exists(paste0("CorMats/", net, "_corMat.csv"))) {
        if (writeOut) {
        }
        if (plotOut) {
            corMat <-
                read.csv(paste0("CorMats/", net, "_corMat.csv"),
                         row.names = 1) %>% as.matrix
            nodes <- getEMSONodes(topoFile)
            corDf <-
                corMat %>% data.frame %>% mutate(Nodes = rownames(.)) %>%
                gather(key = "Nodes1", value = "Correlation",-Nodes) %>%
                mutate(
                    Nodes = factor(Nodes, levels = nodes),
                    Nodes1 = factor(Nodes1, levels = nodes)
                )

            ggplot(corDf, aes(
                x = Nodes,
                y = Nodes1,
                fill = Correlation
            )) + geom_tile() +
                theme_Publication() + scale_fill_gradient2(low = "blue",
                                                           high = "red",
                                                           limits = c(-1, 1)) +
                theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90, vjust = 0.5),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.height = unit(0.5, "in")
                )
            if (!dir.exists("MatrixPlots"))
                dir.create("MatrixPlots")
            ggsave(
                paste0("MatrixPlots/", net, "_corMatPlot.png"),
                width = 7,
                height = 6
            )
            setwd("..")
        }
        if (matOut) {
            corMat <-
                read.csv(paste0("CorMats/", net, "_corMat.csv"),
                         row.names = 1) %>% as.matrix
            return(corMat)
        }
        return()
    }

    net <- str_remove(topoFile, ".topo")
    nodes <-
        readLines(topoFile %>% str_replace(".topo", "_nodes.txt"))
    corMat <-
        read_csv(topoFile %>% str_replace(".topo", "_finFlagFreq.csv"),
                 show_col_types = F) %>%
        filter(flag == 1) %>% select(states, Avg0) %>% drop_na %>%
        mutate(Avg0 = Avg0 * 10000 %>% round) %>% apply(1, function(x) {
            s <- x[1] %>% str_remove_all("'")
            n <- x[2] %>% as.integer
            rep(s, n)
        }) %>% unlist
    if (length(corMat) < 3)
        return(NA)
    corMat <- corMat %>% lapply(function(x) {
        str_split(x, "") %>% unlist %>% as.integer
    }) %>% reduce(rbind.data.frame) %>% set_names(nodes) %>% cor
    colnames(corMat) <- colnames(corMat)
    rownames(corMat) <- rownames(corMat)
    if (plotOut) {
        nodes <- getEMSONodes(topoFile)
        corDf <-
            corMat %>% data.frame %>% mutate(Nodes = rownames(.)) %>%
            gather(key = "Nodes1", value = "Correlation",-Nodes) %>%
            mutate(
                Nodes = factor(Nodes, levels = nodes),
                Nodes1 = factor(Nodes1, levels = nodes)
            )
        textSize <- ifelse(length(nodes) > 28, 0.8, 1)
        if (!is.null(nodeLabelSize))
            textSize <- nodeLabelSize

        ggplot(corDf, aes(
            x = Nodes,
            y = Nodes1,
            fill = Correlation
        )) + geom_tile() +
            theme_Publication() + scale_fill_gradient2(low = "blue",
                                                       high = "red",
                                                       limits = c(-1, 1)) +
            theme(
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.text = element_text(size = rel(textSize)),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.height = unit(0.5, "in")
            ) +
            labs(title = netNameKey[net])
        if (!dir.exists("MatrixPlots"))
            dir.create("MatrixPlots")
        ggsave(
            paste0("MatrixPlots/", net, "_corMatPlot.png"),
            width = 7,
            height = 6
        )
        setwd("..")

    }
    if (writeOut) {
        DirectoryNav("CorMats")
        write.csv(corMat %>% data.frame, paste0(net, "_corMat.csv"))
        setwd("..")
    }

    if (matOut) {
        setwd("..")
        return(corMat)
    }


}
correlationMatBool <- cmpfun(correlationMatBool)

#' Title
#'
#' @param topoFile
#'
#' @return
#' @export
#'
#' @examples
correlationMatrixRAC <- function(topoFile)
{
    wd <- getwd()
    net <- topoFile %>% str_remove(".topo$")
    nodes <- read.delim(paste0(net, ".prs")) %>%
        filter(str_detect(Parameter, "Prod")) %>%
        select(Parameter) %>% unlist %>% str_remove("Prod_of_")
    df <-
        read_delim(paste0(net, "_solution.dat"),
                   delim = "\t",
                   col_names = F) %>%
        set_names(c("ParIndex", "nStates", "Count", nodes)) %>%
        select(-ParIndex,-nStates,-Count)
    corMat <- cor(df)
    nodes <- getEMSONodes(net)
    corDf <-
        corMat %>% data.frame %>% mutate(Nodes = rownames(.)) %>%
        gather(key = "Nodes1", value = "Correlation",-Nodes) %>%
        mutate(
            Nodes = factor(Nodes, levels = nodes),
            Nodes1 = factor(Nodes1, levels = nodes)
        )
    if (plotOut) {
        ggplot(corDf, aes(
            x = Nodes,
            y = Nodes1,
            fill = Correlation
        )) + geom_tile() +
            theme_Publication() + scale_fill_gradient2(low = "blue",
                                                       high = "red",
                                                       limits = c(-1, 1)) +
            theme(
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.height = unit(0.5, "in")
            )
        DirectoryNav("MatrixPlots")
        ggsave(paste0(net, "_corMatPlot.png"),
               width = 7,
               height = 6)
        setwd("..")
    }
    if (writeOut) {
        DirectoryNav("CorMats")
        write.csv(corMat %>% data.frame, paste0(net, "_corMat.csv"))
        setwd("..")
        return()
    }
    if (matOut) {
        setwd("..")
        return(corMat)
    }

}
correlationMatrixRAC <- cmpfun(correlationMatrixRAC)

#' Title
#'
#' @return
#' @export
#'
#' @examples
correlInflDiff <- function() {
    topoFiles <- list.files(".", ".topo$")
    dists <- sapply(topoFiles, function(topoFile)
    {
        # print(topoFile)
        inflFile <-
            paste0("Influence/",
                   str_replace(topoFile, ".topo", "_fullInfl.csv"))
        if (!file.exists(inflFile))
            return(NA)
        inflMat <- read.csv(inflFile, row.names = 1) %>% as.matrix

        corMat <-
            paste0("CorMats/",
                   str_replace(topoFile, ".topo", "_corMat.csv"))
        if (is.na(corMat))
            return(NA)
        inflMat <- inflMat[rownames(corMat), colnames(corMat)]

        sqrt(sum((inflMat - corMat) ^ 2)) / (2 * length(corMat))
    })
    df <-
        data.frame(Network = topoFiles %>% str_remove(".topo"),
                   Dist = dists)
    DirectoryNav("CompiledData")
    write.csv(df, paste0(net, "_correlInflDiff.csv"), row.names = F)
}
correlInflDiff <- cmpfun(correlInflDiff)

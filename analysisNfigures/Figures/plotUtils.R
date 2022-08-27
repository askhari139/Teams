


#' Title
#'
#' @param nets 
#' @param metrics 
#' @param key 
#' @param plotFolder 
#' @param noFlag 
#'
#' @return
#' @export
#'
#' @examples
plotPercentiles <-
    function(nets, metrics, key, plotFolder, noFlag = T) {
        pat <- ifelse(noFlag, "_ALLnoFlag.csv", "_ALL.csv")
        df <- sapply(nets, function(net) {
            setwd(randRaw)
            setwd(net)
            setwd("CompiledData")
            allDat <- read_csv(paste0(net, pat), col_types = cols()) %>%
                mutate(Net = Network) %>% mutate(Network = net)
            gsDat <-
                read_csv(paste0(net, "_TeamStrengths.csv"), col_types = cols()) %>%
                select(Network, Gs) %>% set_names(c("Net", "Gs"))
            allDat <- merge(allDat, gsDat, by = "Net", all = T) %>%
                select(Network, Net, all_of(metrics))
            wtDat <- allDat %>% filter(Net == "wild")
            percentiles <- sapply(metrics, function(x) {
                sum(allDat[[x]] < wtDat[[x]], na.rm = T)
            })
            
            setwd("..")
            percentiles * 100 / nrow(allDat)
            
        })
        if (length(metrics) > 1)
            df <- t(df)
        df <-
            df %>% data.frame %>% set_names(labelshorts[metrics] %>% str_replace_all(" ", "\n")) %>%
            mutate(Network = netNameKey[nets] %>% factor(levels = netNameKey[nets])) %>%
            gather(key = "Metric", value = "Value",-Network) %>%
            mutate(Metric = factor(Metric, levels = labelshorts[metrics] %>% str_replace_all(" ", "\n")))
        ggplot(df, aes(x = Metric, y = Network, fill = Value)) +
            geom_tile() + scale_fill_viridis_c(limits = c(0, 100)) + theme_Publication() +
            theme(
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm")
            ) + labs(x = "", y = "", fill = "Percentile")
        DirectoryNav(plotFolder)
        ggsave(
            paste0(key, "Percentiles.png"),
            width = 2.3 + length(metrics),
            height = 5
        )
    }

#' Title
#'
#' @param nets 
#' @param metric 
#' @param plotFolder 
#' @param noFlag 
#' @param w 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
distribution <- function(nets,
                         metric,
                         plotFolder,
                         noFlag = T,
                         w = NULL,
                         h = NULL) {
    pat <- ifelse(noFlag, "_ALLnoFlag.csv", "_ALL.csv")
    if (metric == "Gs")
        pat <- "_TeamStrengths.csv"
    df <- lapply(nets, function(net) {
        setwd(randRaw)
        setwd(net)
        setwd("CompiledData")
        allDat <- read_csv(paste0(net, pat), col_types = cols()) %>%
            select(Network, all_of(metric)) %>%
            mutate(Net = Network) %>% mutate(Network = net)
        wtDat <- allDat %>% filter(Net == "wild")
        ggplot(allDat, aes_string(x = metric)) +
            geom_histogram(aes(y = ..count.. / sum(..count..))) +
            geom_vline(
                data = wtDat,
                aes_string(xintercept = metric),
                color = "red",
                size = 1.5
            ) +
            labs(x = labelvals[metric],
                 y = "Frequency",
                 title = netNameKey[net]) +
            theme_Publication()
        ggsave(
            paste0(plotFolder, "/", net, "_", metric, ".png"),
            width = 5.5,
            height = 5
        )
        allDat
    }) %>% reduce(rbind.data.frame) %>%
        mutate(Network = netNameKey[Network] %>% factor(levels = netNameKey[nets]))
    wtDf <- df %>% filter(Net == "wild")
    ggplot(df, aes_string(x = metric)) +
        geom_histogram(aes(y = (..count..) / tapply(..count.., ..PANEL.., sum)[..PANEL..])) +
        geom_vline(
            data = wtDf,
            aes_string(xintercept = metric),
            color = "red",
            size = 1.5
        ) +
        facet_wrap( ~ Network, scales = "free", nrow = 1) + theme_Publication() +
        labs(x = labelvals[metric], y = "Frequency")
    h1 <- 5
    w1 <- 3 + 2.5 * length(nets)
    if (!is.null(h))
        h1 <- h
    if (!is.null(w))
        w1 <- w
    
    ggsave(
        paste0(
            plotFolder,
            "/",
            paste0(nets, collapse = "_"),
            "_",
            metric,
            ".png"
        ),
        height = h1,
        width = w1
    )
    
    
}

#' Title
#'
#' @param nets 
#' @param metric 
#' @param plotFolder 
#' @param noFlag 
#' @param h 
#' @param w 
#'
#' @return
#' @export
#'
#' @examples
teamStrengthPlots <-
    function(nets,
             metric,
             plotFolder,
             noFlag = F,
             h = NULL,
             w = NULL) {
        phen <- paste0(metric, "Phen")
        pat <- ifelse(noFlag, "_ALLnoFlag.csv", "_ALL.csv")
        # browser()
        df <- lapply(nets, function(net) {
            # browser()
            setwd(randRaw)
            setwd(net)
            setwd("CompiledData")
            allDat <- read_csv(paste0(net, pat), col_types = cols()) %>%
                select(Network, all_of(c(metric))) %>%
                mutate(Net = Network) %>% mutate(Network = net)
            GsDat <-
                read_csv(paste0(net, "_teamStrengths.csv")) %>% select(Network, Gs) %>%
                set_names(c("Net", "Gs"))
            allDat <-
                merge(allDat, GsDat, by = "Net", all = T) %>% drop_na
            bins <- c(0.2, 0.4, 0.6, 0.8, 1)
            bins <- bins[1:which(bins > max(allDat$Gs))[1]]
            yVal <-
                max(allDat[[metric]]) + 0.1 * (max(allDat[[metric]]) - min(allDat[[metric]]))
            cv <- sapply(bins, function(b) {
                d <- allDat %>% filter(Gs < b & Gs >= (b - 0.2))
                round(sd(d[[metric]], na.rm = T) / mean(d[[metric]], na.rm = T), 2)
            }) %>% data.frame(xVal = bins - 0.1,
                              yVal = yVal,
                              text = .)
            # allDat[[phen]] <- ifelse(allDat[[phen]] == "H", "Hybrid", "Terminal")
            c1 <- cor(allDat$Gs, allDat[[metric]])
            ypos <- ifelse(c1 < -0.1, 0.9, 0.1)
            wtDat <- allDat %>% filter(Net == "wild")
            grob <- correlGrob(allDat, metric, "Gs", yPos = ypos)
            lab <- "WT"
            p <- ggplot(allDat, aes_string(x = "Gs", y = metric)) +
                geom_point(size = 1.5) +
                annotation_custom(grob) +
                geom_text(data = cv, aes(
                    x = xVal,
                    y = yVal,
                    label = text
                ))
            sapply(bins, function(b) {
                p <<-
                    p + geom_vline(
                        xintercept = b,
                        color = "red",
                        linetype = "dashed"
                    )
            })
            p <- p +
                theme_Publication() +
                labs(
                    x = "Team Strength (Ts)",
                    y = labelvals[metric],
                    title = netNameKey[net],
                    color = "Phenotype"
                ) +
                theme(legend.position = "top") +
                geom_label_repel(data = wtDat,
                                 label = lab,
                                 show.legend = F)
            ggsave(
                paste0(plotFolder, "/", net, "_", metric, "_Ts.png"),
                width = 5.5,
                height = 5
            )
            allDat <-
                allDat %>% mutate(Gsbinned = round(Gs, 1) %>% as.character)
            df <-
                allDat %>% select(all_of(metric), Gsbinned) %>% set_names(c("metric", "Gs")) %>%
                filter(is.finite(metric),!is.na(metric),!is.nan(metric))
            aovDat <- aov(metric ~ Gs , data = df)
            pVal <- summary(aovDat)[[1]][["Pr(>F)"]][[1]]
            if (pVal < 0.0001)
                pValKey <- "P-value < 0.0001"
            else if (pVal < 0.001)
                pValKey <- "P-value < 0.001"
            else if (pVal < 0.01)
                pValKey <- "P-value < 0.01"
            else if (pVal < 0.05)
                pValKey <- "P-value < 0.05"
            else
                pValKey <- "n.s."
            allDat$pVal <- pValKey
            # yVal <- min(allDat[[metric]]) + 0.1*(max(allDat[[metric]]) - min(allDat[[metric]]))
            # if (pValKey == "n.s.")
            yVal <-
                max(allDat[[metric]]) + 0.1 * (max(allDat[[metric]]) - min(allDat[[metric]]))
            xVal <- unique(allDat$Gsbinned) %>% sort
            xVal <- xVal[round((length(xVal) + 1) / 2)]
            ggplot(allDat, aes_string(x = "Gsbinned", y = metric)) +
                geom_violin() + theme_Publication() +
                geom_text(
                    aes(x, y, label = caption),
                    data = data.frame(
                        x = xVal,
                        y = yVal,
                        caption = unique(allDat$pVal)
                    ),
                    size = 5.5,
                    fontface = "bold"
                ) +
                labs(x = "Team Strength (Ts)",
                     y = labelvals[metric],
                     title = netNameKey[net])
            ggsave(
                paste0(plotFolder, "/", net, "_", metric, "_TsViolin.png"),
                width = 5.5,
                height = 5
            )
            
            corr <-
                cor.test(allDat[[metric]], allDat[["Gs"]], method = "spearman")
            pVal <- ifelse(corr$p.value < 0.05, "*", "")
            
            
            allDat %>% select(all_of(c("Gs", metric)), Net, Network, Gsbinned, pVal) %>%
                mutate(xVal = xVal, yVal = yVal)
        }) %>% reduce(rbind.data.frame) %>%
            mutate(Network = netNameKey[Network] %>% factor(levels = netNameKey[nets]))
        df$yVal <-
            max(df[[metric]]) + 0.1 * (max(df[[metric]]) - min(df[[metric]]))
        h1 <- 5
        w1 <- 3 + 2.5 * length(nets)
        if (!is.null(h))
            h1 <- h
        if (!is.null(w))
            w1 <- w
        labelDf <- df %>% select(Network, xVal, yVal, pVal) %>% unique
        ggplot(df, aes_string(x = "Gsbinned", y = metric)) +
            geom_violin() + theme_Publication() +
            geom_text(
                aes(x = xVal, y = yVal, label = pVal),
                data = labelDf,
                size = 5.5,
                fontface = "bold"
            ) +
            facet_wrap( ~ Network, nrow = 1, scales = "free_x") +
            labs(x = "Team Strength (Ts)", y = labelvals[metric])
        ggsave(
            paste0(
                plotFolder,
                "/",
                paste0(nets, collapse = "_"),
                "_",
                metric,
                "_Ts.png"
            ),
            height = h1,
            width = w1
        )
    }


#' Title
#'
#' @param nets 
#' @param metrics 
#' @param metric2 
#' @param plotFolder 
#' @param noFlag 
#'
#' @return
#' @export
#'
#' @examples
correlationPlots <-
    function(nets,
             metrics,
             metric2,
             plotFolder,
             noFlag = F) {
        # phen <- paste0(metric, "Phen")
        pat <- ifelse(noFlag, "_ALLnoFlag.csv", "_ALL.csv")
        
        df <- lapply(nets, function(net) {
            print(net)
            setwd(randRaw)
            setwd(net)
            setwd("CompiledData")
            allDat <- read_csv(paste0(net, pat), col_types = cols()) %>%
                # select(Network, all_of(c(metric, phen))) %>%
                mutate(Net = Network) %>% mutate(Network = net)
            GsDat <-
                read_csv(paste0(net, "_teamStrengths.csv")) %>% select(Network, Gs) %>%
                set_names(c("Net", "Gs"))
            allDat <- merge(allDat, GsDat, by = "Net", all = T) %>%
                select(Network, Net, all_of(c(metrics, metric2))) %>% drop_na %>%
                gather(key = "Metric",
                       value = "Value",
                       -all_of(metric2),
                       -Network,
                       -Net) %>%
                mutate(Metric = labelshorts[Metric] %>% str_replace_all(" ", "\n"))
            multiFactorCorrelation(allDat, "Metric", "Value", "Gs", label = F) %>%
                mutate(Network = netNameKey[net])
        }) %>% reduce(rbind.data.frame) %>%
            mutate(
                Label = ifelse(pValue < 0.05, "", "X"),
                Factors = factor(Factors, levels = labelshorts[metrics] %>% str_replace_all(" ", "\n"))
            )
        ggplot(df, aes(x = Factors, y = Network, fill = Correlation)) +
            geom_tile() + geom_text(aes(label = Label)) +
            theme_Publication() +
            scale_x_discrete(expand = c(0, 0)) +
            theme(
                #axis.text.x = element_text(angle = 60, hjust = 1),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm")
            ) +
            scale_fill_gradient2(low = "blue",
                                 high = "red",
                                 limits = c(-1, 1)) +
            labs(
                x = "",
                y = "",
                fill = paste0("Correlation\n", labelshorts[metric2], "\n")
            )
        ggsave(
            paste0(
                plotFolder,
                "/",
                paste0(metrics, collapse = "_"),
                "_",
                metric2,
                ".png"
            ),
            width = 2.3 + length(metrics),
            height = 5
        )
    }

#' Title
#'
#' @param nets 
#' @param metric 
#' @param minMetric 
#' @param maxMetric 
#' @param plotFolder 
#' @param noFlag 
#'
#' @return
#' @export
#'
#' @examples
minMaxViolins <-
    function(nets,
             metric,
             minMetric,
             maxMetric,
             plotFolder,
             noFlag = F) {
        dfAll <- lapply(nets, function(net) {
            # browser()
            setwd(randRaw)
            setwd(net)
            key <- ifelse(noFlag, "noFlag", "")
            randDf <-
                paste0("CompiledData/", net, "_ALL", key, ".csv") %>%
                read_csv(col_types = cols()) %>%
                select(all_of(c(minMetric, maxMetric))) %>%
                gather(key = "Network", value = "Metric") %>%
                mutate(Network = paste0("Rand\n", labelshorts[Network] %>% str_replace_all(" ", "\n")))
            wtDf <- read_csv("wild_finFlagFreq.csv", col_types = cols())
            if (noFlag)
                wtDf <- wtDf %>% filter(flag == 1)
            wtDf <-
                wtDf %>% mutate(Network = "WT") %>% select(Network, all_of(metric)) %>%
                set_names(c("Network", "Metric"))
            df <- rbind.data.frame(wtDf, randDf) %>%
                mutate(Network = factor(
                    Network,
                    levels = c(
                        paste0("Rand ", labelshorts[minMetric]),
                        "WT",
                        paste0("Rand ", labelshorts[maxMetric])
                    ) %>%
                        str_replace_all(" ", "\n")
                ))
            p <-
                ggplot(df, aes(x = Network, y = Metric)) + geom_violin() + theme_Publication() +
                labs(x = "",
                     y = labelshorts[metric],
                     title = netNameKey[net])
            if (metric == "Avg0")
                p <- p + scale_y_log10()
            setwd(plotFolder)
            ggsave(
                paste0(net, "_", labelshorts[metric], "_violin.png"),
                width = 5.5,
                height = 5
            )
            df %>% mutate(Net = netNameKey[net])
        }) %>% reduce(rbind.data.frame)
        p <-
            ggplot(dfAll, aes(x = Network, y = Metric)) + geom_violin() + theme_Publication() +
            labs(x = "", y = labelshorts[metric]) + facet_wrap( ~ Net, nrow = 1)
        if (metric == "Avg0")
            p <- p + scale_y_log10()
        setwd(plotFolder)
        ggsave(
            paste0(paste0(c(
                nets, labelshorts[metric]
            ), collapse = "_"), "_violin.png"),
            width = 3 + 2.5 * length(nets),
            height = 5
        )
        
    }

#' Title
#'
#' @param net 
#'
#' @return
#' @export
#'
#' @examples
singleCoherenceRand <- function(net) {
    setwd(paste0(randRaw, "/", net))
    dfTeams <-
        read_csv(paste0("CompiledData/", net, "_TeamStrengths.csv"))
    setwd("signalCoherence")
    fils <- list.files(".", "_signalCoherence.csv")
    netz <- fils %>% str_remove("_signal.*")
    df <- lapply(fils, function(x) {
        Net <- x %>% str_remove("_signalCoherence.csv")
        d <- read_csv(x)
        groupLabels <-
            readLines("../wild_nodes.txt") %>% str_split(",")
        eNodes <- groupLabels[[1]] %>% length
        mNodes <- groupLabels[[2]] %>% length
        d %>% mutate(
            Core = (Epithelial * eNodes + Mesenchymal * mNodes) / (eNodes + mNodes),
            Network = str_remove(x, "_sign.*")
        )
    }) %>% reduce(rbind.data.frame)
    df <-
        merge(df,
              dfTeams %>% select(Network, Gs),
              by = "Network",
              all.x = T)
    df2 <- df %>%
        group_by(Network, initState, initPhen, signal) %>%
        summarise(
            Epithelial = sum(Freq * Epithelial),
            Mesenchymal = sum(Freq * Mesenchymal),
            Gs = unique(Gs),
            Core = sum(Freq * Core),
            .groups = "drop"
        ) %>%
        group_by(Network, initPhen, signal) %>%
        summarise(
            Epithelial = mean(Epithelial),
            Mesenchymal = mean(Mesenchymal),
            Gs = unique(Gs),
            Core = mean(Core),
            .groups = "drop"
        )
}

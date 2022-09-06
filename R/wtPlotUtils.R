#' Title
#'
#' @param netList 
#' @param plotFolder 
#'
#' @return
#' @export
#'
#' @examples
stabilityPlots <- function(netList, plotFolder)
{
    cohDf <- lapply(netList, function(net) {
        print(net)
        df <-
            read_csv(paste0(randRaw, "/", net, "/wild_finFlagFreq.csv"),
                     col_types = cols()) %>%
            filter(!is.na(Avg0)) %>%
            select(coherence0, Avg0, SD0, frust0) %>%
            mutate(Network = netNameKey[net])
        DirectoryNav(plotFolder)
        ggplot(df, aes(
            x = frust0,
            y = coherence0,
            color = log10(Avg0)
        )) +
            geom_point(size = 1.5) + scale_color_viridis_c() +
            theme_Publication() +
            theme(
                legend.position = c(0.2, 0.4),
                legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                legend.background = element_blank()
            ) +
            labs(x = "Frustration",
                 y = "Coherence",
                 color = "Log10(SSF)")
        ggsave(paste0(net, "_cohVfrust.png"),
               width = 5.5,
               height = 5)
        ggplot(df, aes(
            x = frust0,
            y = Avg0,
            color = coherence0
        )) +
            geom_point(size = 1.5) + scale_color_viridis_c() +
            scale_y_log10() +
            theme_Publication() +
            theme(
                legend.position = c(0.2, 0.4),
                legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                legend.background = element_blank()
            ) +
            labs(x = "Frustration",
                 y = "SSF",
                 color = "Coherence")
        ggsave(paste0(net, "_ssfVfrust.png"),
               width = 5.5,
               height = 5)
        ggplot(df, aes(
            x = Avg0,
            y = coherence0,
            color = frust0
        )) +
            geom_point(size = 1.5) + scale_color_viridis_c() +
            theme_Publication() +
            theme(
                legend.position = c(0.8, 0.4),
                legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                legend.background = element_blank()
            ) +
            labs(x = "SSF",
                 y = "Coherence",
                 color = "Frustration")
        ggsave(paste0(net, "_cohVfreq.png"),
               width = 5.5,
               height = 5)
        
        cr <-
            cor.test(as.numeric(df$coherence0),
                     as.numeric(df$Avg0),
                     method = "spearman")
        df$cohFreq <-
            paste0(
                df$Network,
                ", \u03c1 : ",
                round(cr$estimate, 2),
                ifelse(cr$p.value < 0.05, "*", "")
            )
        cr <-
            cor.test(as.numeric(df$frust0), as.numeric(df$Avg0), method = "spearman")
        df$frustFreq <-
            paste0(
                df$Network,
                ", \u03c1 : ",
                round(cr$estimate, 2),
                ifelse(cr$p.value < 0.05, "*", "")
            )
        cr <-
            cor.test(as.numeric(df$coherence0),
                     as.numeric(df$frust0),
                     method = "spearman")
        df$cohFrust <-
            paste0(
                df$Network,
                ", \u03c1 : ",
                round(cr$estimate, 2),
                ifelse(cr$p.value < 0.05, "*", "")
            )
        df
    }) %>% reduce(rbind.data.frame)
    
    ggplot(cohDf, aes(x = Avg0, y = coherence0, color = cohFreq)) + geom_point(size = 2) +
        # geom_errorbar(aes(ymin = Avg0 - SD0, ymax = Avg0 + SD0)) +
        xlim(min(cohDf$Avg0), 1) +
        scale_x_log10() +
        # geom_smooth(se = F, method = "lm")+
        theme_Publication() + theme(
            legend.position = c(0.3, 0.8),
            legend.direction = "vertical",
            legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = rel(1.1)),
            legend.title = element_text(size = rel(1.1)),
            legend.background = element_blank()
        ) +
        labs(y = "SSF", x = "Coherence", color = "Network")
    setwd(plotFolder)
    ggsave(paste0(paste0(netList, collapse = "_"), "_CohVfreq.png"),
           width = 5.5,
           height = 5)
    ggplot(cohDf, aes(x = frust0, y = coherence0, color = cohFrust)) + geom_point(size = 2) +
        # geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev)) +
        # geom_smooth(se = F, method = "lm")+
        theme_Publication() + theme(
            legend.position = c(0.3, 0.3),
            legend.direction = "vertical",
            legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = rel(1.1)),
            legend.title = element_text(size = rel(1.1)),
            legend.background = element_blank()
        ) +
        labs(x = "Frustration", y = "Coherence", color = "Network")
    setwd(plotFolder)
    ggsave(paste0(paste0(netList, collapse = "_"), "_CohVfrust.png"),
           width = 5.5,
           height = 5)
    ggplot(cohDf, aes(y = frust0, x = Avg0, color = frustFreq)) + geom_point(size = 2) +
        xlim(min(cohDf$Avg0), 1) + ylim(0, 0.5) +
        # geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev)) +
        # geom_smooth(se = F, method = "lm")+
        scale_x_log10() +
        theme_Publication() + theme(
            legend.position = c(0.75, 0.8),
            legend.direction = "vertical",
            legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = rel(1)),
            legend.title = element_text(size = rel(1)),
            # legend.background = element_blank()
        ) +
        labs(y = "Frustration", x = "SSF", color = "Network")
    setwd(plotFolder)
    ggsave(paste0(paste0(netList, collapse = "_"), "_freqVfrust.png"),
           width = 5.5,
           height = 5)
}

#' Title
#'
#' @param net 
#' @param plotFolder 
#' @param w 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
matrixPlot <- function(net,
                       plotFolder,
                       w = NULL,
                       h = NULL)
{
    freqFile <- paste0(randRaw, "/", net, "/wild_finFlagFreq.csv")
    nodes <-
        readLines(paste0(randRaw, "/", net, "/wild_nodes.txt")) %>%
        str_replace_all(regex("\\W+"), "")
    freqData <-
        read_csv(freqFile, col_types = cols()) %>% filter(flag == 1,!is.na(Avg0)) %>%
        select(states, Phenotype, Avg0) %>%
        mutate(Phenotype = factor(Phenotype, levels = rev(c("H", "M", "E")))) %>%
        mutate(states = str_remove_all(states, "'")) %>%
        arrange(Phenotype) %>%
        mutate(states = str_split(states, "") %>% sapply(function(x) {
            paste0(x, collapse = "_")
        })) %>%
        separate(states, nodes, sep = "_")
    nS <- nrow(freqData)
    topoFile <- paste0(randRaw, "/", net, "/wild.topo")
    file.copy(topoFile, paste0(net, ".topo"))
    gr <-
        readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",") %>% unlist
    sig <- nodes[!(nodes %in% gr)]
    nds <- c(gr, sig)
    
    size <- 0.8
    if (length(nodes) > 30)
        size <- 0.6
    gs <-
        readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    topoFile <- paste0(net, ".topo")
    # freqData <- freqData %>%
    #     mutate(Score = freqData %>% select(all_of(gs[[1]])) %>%
    #                mutate_all(as.numeric) %>% rowSums -
    #                freqData %>% select(all_of(gs[[2]])) %>%
    #                mutate_all(as.numeric) %>% rowSums) %>%
    #     mutate(absScore = abs(Score)) %>% arrange(desc(absScore), Score) %>%
    #     mutate(num = paste0(as.character(Phenotype), 1:nrow(.))) %>%
    #     arrange(Avg0*1000 + Score) %>%
    #     mutate(num = factor(num, levels = num))
    freqData <-
        freqData %>% mutate(num = paste0(as.character(Phenotype), 1:nrow(.)))
    # freqGat <- freqData %>%
    #     gather(key = "Nodes", value = "Level", -Phenotype, -num, -Avg0, -Score, -absScore) %>%
    #     mutate(Nodes = factor(Nodes, levels = nds))
    freqGat <- freqData %>%
        gather(key = "Nodes", value = "Level",-Phenotype,-num,-Avg0) %>%
        mutate(Nodes = factor(Nodes, levels = nds)) %>%
        mutate(Level = ifelse(Level == "0", "-1", "1"))
    breaks <- freqData %>%
        mutate(n = 1:nrow(.)) %>%
        split(freqData$Phenotype) %>% sapply(function(x) {
            median(x$n) %>% round
        })
    breaks <- freqData$num[breaks]
    # labels <- c("Terminal", "Hybrid")
    labels <- c("Epithelial", "Mesenchymal", "Hybrid")
    ggplot(freqGat, aes(
        y = Nodes,
        x = reorder(num,-as.numeric(Phenotype)),
        fill = Level
    )) +
        geom_tile(aes(width = 0.1 + Avg0 / max(Avg0))) +
        theme_Publication() +
        theme(
            axis.title.y = element_blank(),
            # axis.text.y = element_text(angle = 90, hjust = 0.5),
            axis.text.y = element_text(
                size = rel(size),
                vjust = 0.5,
                hjust = 1
            ),
            
            legend.position = "right",
            legend.direction = "vertical",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        scale_fill_manual(values = c("white", "black")) +
        scale_x_discrete(breaks = breaks, labels = labels) +
        # scale_y_reverse()+
        labs(x = "", title = netNameKey[net])
    setwd(plotFolder)
    w1 <- 2 + length(nodes) * 0.25
    if (!is.null(w))
        w1 <- w
    h1 <- 4.5 + nS * 0.03
    if (!is.null(h))
        h1 <- h
    ggsave(paste0(net, "_stateFreq.png"),
           width = h1,
           height = w1)
    freqData %>% filter(Phenotype == "H") %>% select(Avg0) %>% unlist %>% sum
}

#' Title
#'
#' @param net 
#' @param plotFolder 
#' @param w 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
signalCoherencePlot <- function(net,
                                plotFolder,
                                w = 10,
                                h = 6) {
    phenKeyRev <- c(Epithelial = "E",
                    Mesenchymal = "M",
                    Hybrid = "H")
    setwd(paste0(randRaw, "/", net, "/signalCoherence"))
    df <- read_csv("wild_signalCoherence.csv", lazy = F)
    df1 <- df %>% gather(key = "Phenotype",
                         value = "Conservation",-initState,
                         -initPhen,
                         -signal,
                         -Freq) %>%
        mutate(initPhen = paste0("Init ", initPhen),
               Phenotype = phenKeyRev[Phenotype])
    ggplot(df1, aes(x = initState, y = signal, fill = Conservation)) +
        geom_tile(height = 0.9, width = 0.9) + facet_grid(Phenotype ~ initPhen, scales = "free") +
        theme_Publication() + scale_fill_viridis_c(limits = c(0, 1)) +
        theme(axis.text.x = element_blank(),
              legend.key.width = unit(0.7, "cm")) +
        labs(
            x = "State",
            y = "Signal",
            title = netNameKey[net],
            fill = "Coherence"
        )
    setwd(plotFolder)
    ggsave(paste0(net, "_signalCoherence.png"),
           width = w,
           height = h)
    
}

coherenceHeatmaps <- function(net)
{#browser()
    setwd(plotFolder)
    if(!dir.exists("coherencePlots"))
        dir.create("coherencePlots")
    setwd("coherencePlots")
    freqFile <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")
    multiCohFile <- paste0(WTcoherence, "/", net, "_coherence.csv")
    nodeCohFile <- paste0(randcompiled, "/", net, "_NodeStateCoherence.csv")
    dfFreq <- read.csv(freqFile)
    dfCoh <- read.csv(multiCohFile)
    dfSingleCoh <- read.csv(nodeCohFile)
    
    
    rownames(dfCoh) <- dfCoh$X
    colnames(dfCoh) <- colnames(dfCoh) %>% str_remove("X") %>% str_remove_all("\\.")
    dfCoh <- dfCoh[, -1] %>% t %>% data.frame %>% mutate(states = rownames(.))
    dfCoh$states <- paste0("'", dfCoh$states, "'")
    df <- merge(dfCoh, dfFreq, by = "states", all.x = T) %>%
        select(-contains("Avg"), -flag, -contains("frust"), -contains("SD"), -Strength, -Partial,
               -Epithelial, -Mesenchymal, -EMTScore, -coherence0)
    df$phenotype <- paste0(df$phenotype, 1:nrow(df))
    
    dfGat <- df %>% gather(key = "nPert", value = "Coherence", -states, -phenotype)
    dfGat$nPert <- dfGat$nPert %>% str_remove("X") %>% as.numeric
    dfGat$nPert <- dfGat$nPert/22
    # dfGat$Phenotype <- paste0(dfGat$phenotype, 1:nrow(dfGat))
    
    ggplot(dfGat, aes(x = nPert, y = phenotype, 
                      fill = Coherence)) + geom_tile() +
        scale_fill_viridis_c(limits = c(0,1)) + theme_Publication() +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_discrete(expand=c(0,0)) + 
        theme(legend.position = "right", legend.direction = "vertical",
              legend.key.height = unit(1, "cm"), 
              legend.key.width = unit(0.6,"cm"),
              legend.title = element_text(hjust = 0.5),
              panel.background = element_blank(), 
              plot.background = element_blank(), panel.grid = element_blank(), 
              axis.text.y = element_blank()) +
        labs(x = "Level of Perturbation", y = "", fill = "Mean\nCoherence")
    ggsave(paste0(net, "_coherenceMultiNode.png"), width = 6.5, height = 5)
    
    
    
    dfCoh <- read.csv(nodeCohFile)
    rownames(dfCoh) <- dfCoh$X
    colnames(dfCoh) <- colnames(dfCoh) %>% str_remove("X") %>% str_remove_all("\\.")
    dfCoh <- dfCoh[, -1] %>% t %>% data.frame %>% mutate(states = rownames(.))
    dfCoh$states <- paste0("'", dfCoh$states, "'")
    df <- merge(dfCoh, dfFreq, by = "states", all.x = T) %>%
        select(-contains("Avg"), -flag, -contains("frust"), -contains("SD"), 
               -Score, -Partial,
               -Epithelial, -Mesenchymal)
    df$phenotype <- paste0(df$phenotype, 1:nrow(df))
    
    dfGat <- df %>% gather(key = "Nodes", value = "Coherence", -states, -phenotype) %>%
        mutate(Nodes = Nodes %>% str_replace_all(regex("\\W+"), ""))
    
    setwd(paste0(randRaw, "/", net))
    l <- groupCalc1(paste0(net, ".topo"))
    nodesC <- l[[1]] %>% unlist
    nodes <- dfGat$Nodes %>% unique
    nodesS <- nodes[!(nodes %in% nodesC)]
    setwd(cwd)
    dfGat$Nodes <- dfGat$Nodes %>% factor(levels = c(nodesC, nodesS))
    ggplot(dfGat, aes(x = Nodes, y = phenotype, 
                      fill = Coherence)) + geom_tile() +
        scale_fill_viridis_c(limits = c(0,1)) + theme_Publication() +
        scale_x_discrete(expand=c(0,0)) + 
        scale_y_discrete(expand=c(0,0)) + 
        theme(legend.position = "right", legend.direction = "vertical",
              legend.key.height = unit(1, "cm"), 
              legend.key.width = unit(0.6,"cm"),
              legend.title = element_text(hjust = 0.5),
              panel.background = element_blank(), 
              plot.background = element_blank(), panel.grid = element_blank(), 
              axis.text.y = element_blank(), 
              axis.text.x = element_text(angle = 90, hjust = 1, size = rel(0.6))) +
        labs(x = "", y = "", fill = "Mean\nCoherence")
    ggsave(paste0(net, "_coherenceSingleNode.png"), width = 6.5, height = 5)
    
    print(net)
}


#' Title
#'
#' @param net
#'
#' @return
#' @export
#'
#' @examples
AllDataFile <- function(net) {
    setwd(randRaw)
    setwd(net)
    metrics <-
        c(
            "minFrust",
            "minFrustPhen",
            "maxFrust",
            "maxFrustPhen",
            "meanFrust",
            "meanNetFrust",
            "minCoh",
            "minCohPhen",
            "maxCoh",
            "maxCohPhen",
            "meanCoh",
            "meanNetCoh",
            "minFreq",
            "minFreqPhen",
            "maxFreq",
            "maxFreqPhen",
            "meanFreq",
            "corFreqFrust",
            "pFreqFrust",
            "corFreqCoh",
            "pFreqCoh",
            "corFreqStren",
            "pFreqStren",
            "corStrenCoh",
            "pStrenCoh",
            "corFrustCoh",
            "pFrustCoh",
            "corFrustStren",
            "pFrustSren",
            "bmSSF",
            "bmCoh",
            "bmFrust",
            "hybridFreq",
            "terminalFreq",
            "nSS"
        )
    topoFiles <- list.files(".", ".topo$")
    df <- lapply(topoFiles, function(topoFile) {
        print(topoFile)
        freqDf <-
            read_csv(
                str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                col_types = cols(),
                lazy = F
            ) %>%
            filter(flag == 1)
        if (nrow(freqDf) < 3)
            return(rep(NA, 35))
        frust <- freqDf$frust0
        freq <- freqDf$Avg0
        coh <- freqDf$coherence0
        stren <- freqDf$Strength
        if (is.null(coh)) {
            print(topoFile)
            freqDf <-
                read_csv(
                    str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                    col_types = cols(),
                    lazy = F
                )
            d <- coherence(topoFile, write = F)
            if (is.na(d)) {
                freqDf$coherence0 <- NA
            }
            else {
                d <-
                    d[d$init == d$fin, ] %>% mutate(states = init, coherence0 = Freq) %>%
                    select(states, coherence0)
                coherenceVec <- d$coherence0
                names(coherenceVec) <- d$states

                freqDf$coherence0 <- coherenceVec[freqDf$states]
            }

            write_csv(
                freqDf,
                str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                quote = "none",
                na = ""
            )
            coh <- freqDf %>% filter(flag == 1) %>%
                select(coherence0) %>% unlist
        }
        frustration <-
            c(
                frust %>% min(na.rm = T),
                freqDf$Phenotype[which.min(frust)],
                frust %>% max(na.rm = T),
                freqDf$Phenotype[which.max(frust)],
                frust %>% mean(na.rm = T),
                sum(frust * freq, na.rm = T)
            )
        coherence <-
            c(
                coh %>% min(na.rm = T),
                freqDf$Phenotype[which.min(coh)],
                coh %>% max(na.rm = T),
                freqDf$Phenotype[which.max(coh)],
                coh %>% mean(na.rm = T),
                sum(coh * freq, na.rm = T)
            )
        frequency <-
            c(
                freq %>% min(na.rm = T),
                freqDf$Phenotype[which.min(freq)],
                freq %>% max(na.rm = T),
                freqDf$Phenotype[which.max(freq)],
                freq %>% mean(na.rm = T)
            )
        cohFreq <-
            cohFrust <- cohStren <- list(estimate = NA, p.value = NA)
        if (!is.na(coh)) {
            cohFreq <- cor.test(coh, log10(freq), method = "spearman")
            cohFrust <- cor.test(coh, frust, method = "spearman")
            cohStren <- cor.test(coh, stren, method = "spearman")
        }

        freqFrust <-
            cor.test(log10(freq), frust, method = "spearman")
        freqStren <-
            cor.test(log10(freq), stren, method = "spearman")
        frustStren <- cor.test(frust, stren, method = "spearman")

        cors <- c(
            freqFrust$estimate,
            freqFrust$p.value,
            cohFreq$estimate,
            cohFreq$p.value,
            freqStren$estimate,
            freqStren$p.value,
            cohStren$estimate,
            cohStren$p.value,
            cohFrust$estimate,
            cohFrust$p.value,
            frustStren$estimate,
            frustStren$p.value
        )
        bimodalities <- c(
            bimodality_coefficient(log10(freq)),
            bimodality_coefficient(frust),
            bimodality_coefficient(coh)
        )
        hybridFreq <-
            freqDf %>% filter(Phenotype == "H") %>% select(Avg0) %>%
            unlist %>% sum
        terminalFreq <-
            freqDf %>% filter(Phenotype %in% c("E", "M")) %>% select(Avg0) %>%
            unlist %>% sum
        c(
            frustration,
            coherence,
            frequency,
            cors,
            bimodalities,
            hybridFreq,
            terminalFreq,
            length(freq)
        )
    }) %>% reduce(rbind.data.frame) %>% set_names(metrics) %>%
        mutate(Network = topoFiles %>% str_remove(".topo"))
    DirectoryNav("CompiledData")
    write_csv(df, paste0(net, "_ALL.csv"))
    setwd("..")
}

#' Title
#'
#' @param net
#'
#' @return
#' @export
#'
#' @examples
AllDataFileNoFlag <- function(net) {
    setwd(randRaw)
    setwd(net)
    metrics <-
        c(
            "minFrust",
            "minFrustPhen",
            "maxFrust",
            "maxFrustPhen",
            "meanFrust",
            "meanNetFrust",
            "minCoh",
            "minCohPhen",
            "maxCoh",
            "maxCohPhen",
            "meanCoh",
            "meanNetCoh",
            "minFreq",
            "minFreqPhen",
            "maxFreq",
            "maxFreqPhen",
            "meanFreq",
            "corFreqFrust",
            "pFreqFrust",
            "corFreqCoh",
            "pFreqCoh",
            "corFreqStren",
            "pFreqStren",
            "corStrenCoh",
            "pStrenCoh",
            "corFrustCoh",
            "pFrustCoh",
            "corFrustStren",
            "pFrustSren",
            "bmSSF",
            "bmCoh",
            "bmFrust",
            "hybridFreq",
            "terminalFreq",
            "nSS"
        )
    topoFiles <- list.files(".", ".topo$")
    df <- lapply(topoFiles, function(topoFile) {
        # print(topoFile)
        freqDf <-
            read_csv(
                str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                col_types = cols(),
                lazy = F
            )
        if (nrow(freqDf) < 3)
            return(rep(NA, 35))
        frust <- freqDf$frust0
        freq <- freqDf$Avg0
        coh <- freqDf$coherence0
        stren <- freqDf$Strength
        if (is.null(coh)) {
            print(topoFile)
            freqDf <-
                read_csv(
                    str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                    col_types = cols(),
                    lazy = F
                )
            d <- coherence(topoFile, write = F)
            if (is.na(d)) {
                freqDf$coherence0 <- NA
            }
            else {
                d <-
                    d[d$init == d$fin, ] %>% mutate(states = init, coherence0 = Freq) %>%
                    select(states, coherence0)
                coherenceVec <- d$coherence0
                names(coherenceVec) <- d$states

                freqDf$coherence0 <- coherenceVec[freqDf$states]
            }

            write_csv(
                freqDf,
                str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                quote = "none",
                na = ""
            )
            coh <- freqDf %>% filter(flag == 1) %>%
                select(coherence0) %>% unlist
        }
        frustration <-
            c(
                frust %>% min(na.rm = T),
                freqDf$Phenotype[which.min(frust)],
                frust %>% max(na.rm = T),
                freqDf$Phenotype[which.max(frust)],
                frust %>% mean(na.rm = T),
                sum(frust * freq, na.rm = T)
            )
        coherence <-
            c(
                coh %>% min(na.rm = T),
                freqDf$Phenotype[which.min(coh)],
                coh %>% max(na.rm = T),
                freqDf$Phenotype[which.max(coh)],
                coh %>% mean(na.rm = T),
                sum(coh * freq, na.rm = T)
            )
        frequency <-
            c(
                freq %>% min(na.rm = T),
                freqDf$Phenotype[which.min(freq)],
                freq %>% max(na.rm = T),
                freqDf$Phenotype[which.max(freq)],
                freq %>% mean(na.rm = T)
            )
        cohFreq <-
            cohFrust <- cohStren <- list(estimate = NA, p.value = NA)
        if (!is.na(coh)) {
            cohFreq <- cor.test(coh, log10(freq), method = "spearman")
            cohFrust <- cor.test(coh, frust, method = "spearman")
            cohStren <- cor.test(coh, stren, method = "spearman")
        }

        freqFrust <-
            cor.test(log10(freq), frust, method = "spearman")
        freqStren <-
            cor.test(log10(freq), stren, method = "spearman")
        frustStren <- cor.test(frust, stren, method = "spearman")

        cors <- c(
            freqFrust$estimate,
            freqFrust$p.value,
            cohFreq$estimate,
            cohFreq$p.value,
            freqStren$estimate,
            freqStren$p.value,
            cohStren$estimate,
            cohStren$p.value,
            cohFrust$estimate,
            cohFrust$p.value,
            frustStren$estimate,
            frustStren$p.value
        )
        bimodalities <- c(
            bimodality_coefficient(log10(freq)),
            bimodality_coefficient(frust),
            bimodality_coefficient(coh)
        )
        hybridFreq <-
            freqDf %>% filter(Phenotype == "H") %>% select(Avg0) %>%
            unlist %>% sum
        terminalFreq <-
            freqDf %>% filter(Phenotype %in% c("E", "M")) %>% select(Avg0) %>%
            unlist %>% sum
        c(
            frustration,
            coherence,
            frequency,
            cors,
            bimodalities,
            hybridFreq,
            terminalFreq,
            length(freq)
        )
    }) %>% reduce(rbind.data.frame) %>% set_names(metrics) %>%
        mutate(Network = topoFiles %>% str_remove(".topo"))
    DirectoryNav("CompiledData")
    write_csv(df, paste0(net, "_ALLnoFlag.csv"))
    setwd("..")
}

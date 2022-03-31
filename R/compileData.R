AllDataFile <- function(net) {
    metrics <- c("minFrust", "maxFrust", "meanFrust", "meanNetFrust",
                 "minCoh", "maxCoh", "meanCoh", "meanNetCoh",
                 "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh",
                 "corFrustCoh", "pFrustCoh", "bmSSF", "bmCoh", "bmFrust")
    topoFiles <- list.files(".", ".topo$")
    df <- sapply(topoFiles, function(topoFile) {
        freqDf <- read_csv(str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                           col_types = cols(), lazy = F) %>% filter(flag == 1)
        if (nrow(freqDf) < 3)
            return(rep(NA, 17))
        frust <- freqDf$frust0
        freq <- freqDf$Avg0
        coh <- freqDf$coherence0
        frustration <- c(frust %>% min(na.rm = T), frust %>% max(na.rm = T),
                        frust %>% mean(na.rm = T), sum(frust*freq,na.rm = T))
        coherence <- c(coh %>% min(na.rm = T), coh %>% max(na.rm = T),
                       coh %>% mean(na.rm = T), sum(coh*freq,na.rm = T))
        cohFreq <- cor.test(coh, freq, method = "spearman")
        cohFrust <- cor.test(coh, frust, method = "spearman")
        freqFrust <- cor.test(freq, frust, method = "spearman")
        cors <- c(freqFrust$estimate, freqFrust$p.value,
                  cohFreq$estimate, cohFreq$p.value,
                  cohFrust$estimate, cohFrust$p.value)
        bimodalities <- c(bimodality_coefficient(freq),
                          bimodality_coefficient(frust),
                          bimodality_coefficient(coh))
        c(frustration, coherence, cors, bimodalities)
    }) %>% t %>% data.frame %>% set_names(metrics) %>%
        mutate(Network = topoFiles %>% str_remove(".topo"))
    DirectoryNav("CompiledData")
    write_csv(df, paste0(net, "_ALL.csv"))
    setwd("..")
}

AllDataFileNoFlag <- function(net) {
    metrics <- c("minFrust", "maxFrust", "meanFrust", "meanNetFrust",
                 "minCoh", "maxCoh", "meanCoh", "meanNetCoh",
                 "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh",
                 "corFrustCoh", "pFrustCoh", "bmSSF", "bmCoh", "bmFrust")
    topoFiles <- list.files(".", ".topo$")
    df <- sapply(topoFiles, function(topoFile) {
        freqDf <- read_csv(str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
                           col_types = cols(), lazy = F)
        if (nrow(freqDf) < 3)
            return(rep(NA, 17))
        frust <- freqDf$frust0
        freq <- freqDf$Avg0
        coh <- freqDf$coherence0
        frustration <- c(frust %>% min(na.rm = T), frust %>% max(na.rm = T),
                         frust %>% mean(na.rm = T), sum(frust*freq,na.rm = T))
        coherence <- c(coh %>% min(na.rm = T), coh %>% max(na.rm = T),
                       coh %>% mean(na.rm = T),sum(coh*freq, na.rm = T))
        cohFreq <- cor.test(coh, freq, method = "spearman")
        cohFrust <- cor.test(coh, frust, method = "spearman")
        freqFrust <- cor.test(freq, frust, method = "spearman")
        cors <- c(freqFrust$estimate, freqFrust$p.value,
                  cohFreq$estimate, cohFreq$p.value,
                  cohFrust$estimate, cohFrust$p.value)
        bimodalities <- c(bimodality_coefficient(freq),
                          bimodality_coefficient(frust),
                          bimodality_coefficient(coh))
        c(frustration, coherence, cors, bimodalities)
    }) %>% t %>% data.frame %>% set_names(metrics) %>%
        mutate(Network = topoFiles %>% str_remove(".topo"))
    DirectoryNav("CompiledData")
    write_csv(df, paste0(net, "_ALLNoflag.csv"))
    setwd("..")
}

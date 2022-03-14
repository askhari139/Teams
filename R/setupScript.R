installed <- library()$results[, 1]
packages <- c("tidyr", "readr", "magrittr", "ggplot2", "compiler","dplyr", "JuliaCall", "xfun")
packages <- packages[!(packages %in% installed)]

library(tidyverse)
library(compiler)
library(ggrepel)
library(grid)
library(ggcorrplot)
options(stringsAsFactors = F, lazy = F)

mainFolder <- readline(prompt="Enter the path of the folder to save the data: ")
topoFolder <- readline(prompt="Enter the path of the folder that has the topofiles: ")
numThreads <- readline(prompt = "Enter the number of threads available for simulation: ") %>% as.integer
if(is.na(numThreads))
    numThreads <- 1

Sys.setenv(JULIA_NUM_THREADS = as.character(numThreads))

### Input files----
topoFileFolder <- paste0(mainFolder, "/TopoFiles")

### Simulation Package----
simPackage <- paste0(mainFolder, "/BooleanSim")

### Analysis Data folders-----
results <- paste0(mainFolder, "/Final_Results")
randcompiled <- paste0(results, "/Compiled_data")
randRaw <- paste0(results,"/Raw_Data")
WTcoherence <- paste0(results, "/WT_NPerts")
QC <- paste0(results, "/QC")
RACIPE <- paste0(results,"/RACIPE")
RACIPE_WT <- paste0(RACIPE, "/WildType")
RACIPE_QC <- paste0(RACIPE, "/QC")
gsCausation <- paste0(mainFolder, "/Final_Results/GsMultiPert")

### Figure folder----
cwd <- paste0(mainFolder, "/Figures")

### Misc folders----
edgeDel <- paste0(mainFolder, "/Final_Results/edge_deletions")
RACIPE_OE <- paste0(RACIPE, "/OverExpression")
RACIPE_DE <- paste0(RACIPE, "/DownExpression")
NodeDel <- paste0(results, "/NodePert")

folderList <- c(mainFolder, topoFileFolder, simPackage, results, randcompiled, randRaw, WTcoherence,
                QC, RACIPE, RACIPE_WT, RACIPE_QC,
                cwd, edgeDel, RACIPE_OE, RACIPE_DE, NodeDel, gsCausation)
sapply(folderList, function(fol) {
    if (!dir.exists(fol))
        dir.create(fol)
})

setwd(topoFolder)
topoFiles <- list.files(".", ".topo$")
sapply(topoFiles, function(topoFile) {
    file.copy(topoFile, topoFileFolder)
})

netList <- topoFiles %>% str_remove(".topo")

# netList <- c("EMT_RACIPE", "EMT_RACIPE2", "SIL", "SIL2",
# "melanoma", "SCLC", "EMT_MET_reduced", "drosophila")
EMPNets <- netList

labelKeys <- c("maxFrust", "minFrust", "meanFrust", "maxNetFrust", "minNetFrust",
               "meanNetFrust","maxCoh", "minCoh", "meanCoh", "maxNetCoh", "minNetCoh",
               "meanNetCoh", "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh",
               "corFrustCoh", "pFrustCoh"   )
labelvals <- c("Maximum Frustration", "Minimum Frustration", "Mean Frustration",
               "Max Net Frustration", "Min Net Frustration", "Mean Net Frustration",
               "Maximum Coherence", "Minimum Coherence", "Mean Coherence",
               "Max Net Coherence", "Min Net Coherence", "Mean Net Coherence",
               "\u03c1 Frequency-frustration", "pVal Frequency-frustration",
               "\u03c1 Frequency-coherence", "pVal Frequency-coherence",
               "\u03c1 Coherence-frustration", "pVal Coherence-frustration")
names(labelvals) <- labelKeys

labelshorts <- c("Max Frust", "Min Frust", "Mean Frust",
                 "Max Net Frust", "Min Net Frust", "Mean Net Frust",
                 "Max Coh", "Min Coh", "Mean Coh",
                 "Max Net Coh", "Min Net Coh", "Mean Net Coh",
                 "\u03c1 Freq-Frust", "pVal Freq-Frust",
                 "\u03c1 Freq-Coh", "pVal Freq-Coh",
                 "\u03c1 Coh-Frust", "pVal Coh-Frust")
names(labelshorts) <- labelKeys

setwd(topoFileFolder)
netNameKey <- sapply(topoFiles, function(x){
    df <- read.delim(x, sep = "") %>%
        mutate_if(is.character, str_replace_all, pattern = regex("\\W+"), replace = "")
    nodes <- unique(c(df$Source, df$Target)) %>% length
    edges <- nrow(df)
    write_delim(df, x, delim = " ", quote = "none")
    paste0(nodes, "N ", edges, "E")
}) %>% set_names(netList)

netNameKeyRAC <- c(EMT_RACIPE = "22N 82E", EMT_RACIPE2 = "26N 100E",
                   silviera = "18N 33E", silviera2 = "20N 40E", EMT_MET_reduced = "57N 113E")
phenKey <- c(E = "Epithelial", H = "Hybrid", M = "Mesenchymal")


# source(paste0(cwd,"/plot_theme.R"))
# source(paste0(cwd,"/utils.R"))
# source(paste0(mainFolder, "\\Final_Results\\codes\\inflMat.R"))
# source(paste0(mainFolder, "\\Final_Results\\codes\\groupPlotter.R"))

BmodelSetup <- function() {
    setwd(simPackage)
    download.file("https://github.com/askhari139/Boolean.jl/archive/refs/heads/main.zip",
                  destfile = "Bmodel.zip")
    zipFile <- "Bmodel.zip"
    unzip(zipFile, junkpaths = T)
    file.remove(zipFile)
    script <- readLines("script.jl")
    script[1] <- paste0("include(\"", simPackage, "/bmodel.jl\")")
    writeLines(script, "script.jl")
    julia_source("dependencyInstaller.jl")
}
BmodelSetup()

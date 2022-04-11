installed <- library()$results[, 1]
packages <- c("tidyr", "readr", "magrittr", "ggplot2", "compiler", "dplyr",
              "JuliaCall", "xfun", "stringr", "future", "future.apply", "purrr",
              "devtools")
toInstall <- packages[!(packages %in% installed)]
sapply(toInstall, install.packages, repos = "https://cloud.r-project.org/")

if(!require(miscFuncs))
  devtools::install_github("askhari139/miscFuncs")
packages <- c(packages, "miscFuncs")

sapply(packages, library, character.only = T)
options(stringsAsFactors = F, lazy = F)


# Sys.setenv(JULIA_NUM_THREADS = as.character(numThreads))


SetupFunc <- function(mainFolder, topoFolder, numThreads = 1) {
  ### Input files----
  topoFileFolder <<- paste0(mainFolder, "/TopoFiles")

  ### Simulation Package----
  simPackage <<- paste0(mainFolder, "/BooleanSim")

  ### Analysis Data folders-----
  results <<- paste0(mainFolder, "/Final_Results")
  randcompiled <<- paste0(results, "/Compiled_data")
  randRaw <<- paste0(results, "/Raw_Data")
  WTcoherence <<- paste0(results, "/WT_NPerts")
  QC <<- paste0(results, "/QC")
  RACIPE <<- paste0(results, "/RACIPE")
  RACIPE_WT <<- paste0(RACIPE, "/WildType")
  RACIPE_QC <<- paste0(RACIPE, "/QC")
  gsCausation <<- paste0(mainFolder, "/Final_Results/GsMultiPert")

  ### Figure folder----
  cwd <<- paste0(mainFolder, "/Figures")

  ### Misc folders----
  edgeDel <<- paste0(mainFolder, "/Final_Results/edge_deletions")
  RACIPE_OE <<- paste0(RACIPE, "/OverExpression")
  RACIPE_DE <<- paste0(RACIPE, "/DownExpression")
  NodeDel <<- paste0(results, "/NodePert")

  folderList <<- c(
    mainFolder, topoFileFolder, simPackage, results, randcompiled, randRaw, WTcoherence,
    QC, RACIPE, RACIPE_WT, RACIPE_QC,
    cwd, edgeDel, RACIPE_OE, RACIPE_DE, NodeDel, gsCausation
  )
  sapply(folderList, function(fol) {
    if (!dir.exists(fol)) {
      dir.create(fol)
    }
  })

  setwd(topoFolder)
  topoFiles <<- list.files(".", ".topo$")
  sapply(topoFiles, function(topoFile) {
    file.copy(topoFile, topoFileFolder)
  })

  netList <<- topoFiles %>% str_remove(".topo")

  # netList <<- c("EMT_RACIPE", "EMT_RACIPE2", "SIL", "SIL2",
  # "melanoma", "SCLC", "EMT_MET_reduced", "drosophila")
  EMPNets <<- netList

  labelKeys <<- c(
    "maxFrust", "minFrust", "meanFrust", "maxNetFrust", "minNetFrust",
    "meanNetFrust", "maxCoh", "minCoh", "meanCoh", "maxNetCoh", "minNetCoh",
    "meanNetCoh", "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh",
    "corFrustCoh", "pFrustCoh", "Gs", "bmSSF", "bmFrust", "bmCoh"
  )
  labelvals <<- c(
    "Maximum Frustration", "Minimum Frustration", "Mean Frustration",
    "Max Net Frustration", "Min Net Frustration", "Mean Net Frustration",
    "Maximum Coherence", "Minimum Coherence", "Mean Coherence",
    "Max Net Coherence", "Min Net Coherence", "Mean Net Coherence",
    "\u03c1 Frequency-frustration", "pVal Frequency-frustration",
    "\u03c1 Frequency-coherence", "pVal Frequency-coherence",
    "\u03c1 Coherence-frustration", "pVal Coherence-frustration",
    "Mean Team Strength", "SSF Bimodality Coefficient",
    "Frustration Bimodality Coefficient", "Coherence Bimodality Coefficient"
  )
  names(labelvals) <<- labelKeys

  labelshorts <<- c(
    "Max Frust", "Min Frust", "Mean Frust",
    "Max Net Frust", "Min Net Frust", "Mean Net Frust",
    "Max Coh", "Min Coh", "Mean Coh",
    "Max Net Coh", "Min Net Coh", "Mean Net Coh",
    "\u03c1 Freq-Frust", "pVal Freq-Frust",
    "\u03c1 Freq-Coh", "pVal Freq-Coh",
    "\u03c1 Coh-Frust", "pVal Coh-Frust", "Team Strength", "SSF Bimodality",
    "Frust Bimodality", "Coh Bimodality"
  )
  names(labelshorts) <<- labelKeys

  setwd(topoFileFolder)
  netNameKey <<- sapply(topoFiles, function(x) {
    df <- read.delim(x, sep = "", row.names = NULL) %>%
      mutate_if(is.character, str_replace_all, pattern = regex("\\W+"), replace = "") %>%
      mutate_if(is.character, tolower)
    nodes <- unique(c(df$Source, df$Target)) %>% length()
    edges <- nrow(df)
    write_delim(df, x, delim = " ", quote = "none")
    paste0(nodes, "N ", edges, "E")
  }) %>% set_names(netList)

  netNameKeyRAC <<- c(
    EMT_RACIPE = "22N 82E", EMT_RACIPE2 = "26N 100E",
    silviera = "18N 33E", silviera2 = "20N 40E", EMT_MET_reduced = "57N 113E"
  )
  phenKey <<- c(E = "Epithelial", H = "Hybrid", M = "Mesenchymal")
  mainFolder <<- mainFolder
  numThreads <<- numThreads
  BmodelSetup()
}




BmodelSetup <- function() {
  setwd(simPackage)
  download.file("https://github.com/askhari139/Boolean.jl/archive/refs/heads/main.zip",
    destfile = "Bmodel.zip"
  )
  zipFile <- "Bmodel.zip"
  unzip(zipFile, junkpaths = T)
  file.remove(zipFile)
  script <- readLines("script.jl")
  script[1] <- paste0("include(\"", simPackage, "/bmodel.jl\")")
  writeLines(script, "script.jl")

  script <- readLines("scriptWindows.jl")
  script[1] <- paste0("include(\"", simPackage, "/bmodel.jl\")")
  writeLines(script, "scriptWindows.jl")
  julia_source("dependencyInstaller.jl")
}

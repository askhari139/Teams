installed <- library()$results[, 1]
packages <-
    c(
        "tidyr",
        "readr",
        "magrittr",
        "ggplot2",
        "compiler",
        "dplyr",
        "JuliaCall",
        "xfun",
        "stringr",
        "future",
        "future.apply",
        "purrr",
        "devtools"
    )
toInstall <- packages[!(packages %in% installed)]
sapply(toInstall, install.packages, repos = "https://cloud.r-project.org/")

if (!require(funcsKishore))
    devtools::install_github("askhari139/funcsKishore")
packages <- c(packages, "funcsKishore")

sapply(packages, library, character.only = T)
options(stringsAsFactors = F, lazy = F)


# Sys.setenv(JULIA_NUM_THREADS = as.character(numThreads))


#' Title
#'
#' @param mainFolder
#' @param topoFolder
#' @param numThreads
#' @param julia
#'
#' @return
#' @export
#'
#' @examples
SetupFunc <-
    function(mainFolder,
             topoFolder,
             numThreads = 1,
             julia = F) {
        ### Input files----
        topoFileFolder <<- paste0(mainFolder, "/TopoFiles")

        ### Simulation Package----
        simPackage <<- paste0(mainFolder, "/BooleanSim")

        ### Analysis Data folders-----
        results <<- paste0(mainFolder, "/Final_Results")
        randRaw <<- paste0(results, "/Raw_Data")
        WTcoherence <<- paste0(results, "/WT_NPerts")
        RACIPE <<- paste0(results, "/RACIPE")
        RACIPE_WT <<- paste0(RACIPE, "/WildType")
        gsCausation <<- paste0(mainFolder, "/Final_Results/GsMultiPert")

        ### Figure folder----
        cwd <<- paste0(mainFolder, "/Figures")

        ### Misc folders----
        # edgeDel <<- paste0(mainFolder, "/Final_Results/edge_deletions")
        # RACIPE_OE <<- paste0(RACIPE, "/OverExpression")
        # RACIPE_DE <<- paste0(RACIPE, "/DownExpression")
        # NodeDel <<- paste0(results, "/NodePert")

        folderList <<- c(
            mainFolder,
            topoFileFolder,
            simPackage,
            results,
            randRaw,
            RACIPE,
            RACIPE_WT,
            cwd,
            # edgeDel,
            # RACIPE_OE,
            # RACIPE_DE,
            # NodeDel,
            gsCausation
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
        EMPNets <<- c("SIL",
                      "SIL2",
                      "EMT_RACIPE",
                      "EMT_RACIPE2",
                      "EMT_MET_reduced")

        extras <<- c("sclcnetwork", "melanoma")

        labelKeys <<- c(
            "Avg0",
            "frust0",
            "coherence0",
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
            "corFrustCoh",
            "pFrustCoh",
            "bmSSF",
            "bmCoh",
            "bmFrust",
            "hybridFreq",
            "terminalFreq",
            "Gs",
            "nSS"
        )
        labelvals <<- c(
            "SSF",
            "Frustration",
            "Coherence",
            "Minimum Frustration",
            "Minimum Frustration Phenotype",
            "Maximum Frustration",
            "Maximum Frustration Phenotype",
            "Mean Frustration",
            "Mean Net Frustration",
            "Minimum Coherence",
            "Minimum Coherence Phenotype",
            "Maximum Coherence",
            "Maximum Coherence Phenotype",
            "Mean Coherence",
            "Mean Net Coherence",
            "Minimum SSF",
            "Minimum SSF Phenotype",
            "Maximum SSF",
            "Maximum SSF Phenotype",
            "Mean SSF",
            "\u03c1 SSF-frustration",
            "pVal SSF-frustration",
            "\u03c1 SSF-coherence",
            "pVal SSF-coherence",
            "\u03c1 Coherence-frustration",
            "pVal Coherence-frustration",
            "SSF Bimodality Coefficient",
            "Coherence Bimodality Coefficient",
            "Frustration Bimodality Coefficient",
            "Hybrid",
            "Terminal",
            "Team Strength (Ts)",
            "# Steady States"
        )
        names(labelvals) <<- labelKeys

        labelshorts <<- c(
            "SSF",
            "Frustration",
            "Coherence",
            "Min Frust",
            "Min Frust Phenotype",
            "Max Frust",
            "Max Frust Phenotype",
            "Mean Frust",
            "Mean Net Frustration",
            "Min Coherence",
            "Min Coherence Phenotype",
            "Max Coherence",
            "Max Coherence Phenotype",
            "Mean Coherence",
            "Mean Net Coherence",
            "Min SSF",
            "Min SSF Phenotype",
            "Max SSF",
            "Max SSF Phenotype",
            "Mean SSF",
            "\u03c1 SSF Frust",
            "pVal SSF-Frust",
            "\u03c1 SSF Coh",
            "pVal SSF Coh",
            "\u03c1 Coh Frust",
            "pVal Coh-Frust",
            "B.C. SSF",
            "B.C. Coherence",
            "B.C. Frust",
            "Hybrid",
            "Terminal",
            "Ts",
            "Num SS"
        )
        names(labelshorts) <<- labelKeys

        setwd(topoFileFolder)
        netNameKey <<- sapply(topoFiles, function(x) {
            df <- read.delim(x, sep = "", row.names = NULL) %>%
                mutate_if(
                    is.character,
                    str_replace_all,
                    pattern = regex("\\W+"),
                    replace = ""
                )
            nodes <- unique(c(df$Source, df$Target)) %>% length()
            edges <- nrow(df)
            write_delim(df, x, delim = " ", quote = "none")
            paste0(nodes, "N ", edges, "E")
        }) %>% set_names(netList)

        netNameKeyRAC <<- c(
            EMT_RACIPE = "22N 82E",
            EMT_RACIPE2 = "26N 100E",
            silviera = "18N 33E",
            silviera2 = "20N 40E",
            EMT_MET_reduced = "57N 113E"
        )
        phenKey <<- c(E = "Epithelial", H = "Hybrid", M = "Mesenchymal")
        mainFolder <<- mainFolder
        numThreads <<- numThreads
        if (julia)
            BmodelSetup()
    }




BmodelSetup <- function() {
    setwd(simPackage)
    download.file(
        "https://github.com/askhari139/Boolean.jl/archive/refs/heads/main.zip",
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

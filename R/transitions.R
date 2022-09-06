## sigmoidal Fits

sigmFit <- function(x)
{
    df <- read_csv(x, show_col_types = F)
    if (ncol(df) < 3)
        return()
    df <- df %>% select(init, initPhen, nNode, Hamming, Freq) %>%
        mutate(Fraction = nNode / max(nNode)) %>%
        group_by(init, Fraction) %>%
        summarise(
            Hamming = sum(Hamming * Freq),
            initPhen = unique(initPhen),
            .groups = "drop"
        )
    inits <- unique(df$init)

    fits <- sapply(inits, function(x) {
        d <- df %>% filter(init == x)
        xDat <- log(0.5 / d$Fraction)
        yDat <- log((1 + 0.002) / (d$Hamming + 0.001) - 1)
        fit <- lm(yDat ~ xDat)
        c(fit$coefficients[2],
          fit$coefficients[1],
          unique(d$initPhen))
    }) %>% t %>% data.frame %>% set_names(c("Cooperativity", "Intercept", "Phenotype")) %>%
        mutate(Net = x %>% str_remove("_allNodeCoherence.*"),
               State = inits)
    return(fits)
}
#' Title
#'
#' @param net
#'
#' @return
#' @export
#'
#' @examples
sigmDat <- function(net)
{
    setwd(paste0(randRaw, "/", net, "/PhenotypicTransition"))
    # fitsWT <- paste0("wild_allNodeCoherence_nPert100_nIter10_reps1.csv") %>%
    #     sigmFit %>% mutate(Net = net)
    filz <- list.files(".", "allNodeCoherence")
    plan(multisession, workers = 8)
    fittedDat <-
        future_lapply(filz, sigmFit) %>% reduce(rbind.data.frame)
    future:::ClusterRegistry("stop")
    # setwd("..")
    setwd("../CompiledData")
    df <- fittedDat
    write_csv(df, paste0(net, "_sigmoidalFits.csv"))
    print(net)
}


emma.delta.ML.dLL.wo.Z <- function (logdelta, lambda, etas, xi)
{
    n <- length(xi)
    delta <- exp(logdelta)
    etasq <- etas * etas
    ldelta <- lambda + delta
    return(0.5 * (n * sum(etasq/(ldelta * ldelta))/sum(etasq/ldelta) -
        sum(1/(xi + delta))))
}








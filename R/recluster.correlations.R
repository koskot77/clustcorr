recluster.correlations <- function(m,cutoff){
    .C("reclusterCorrelations",as.numeric(cutoff),clust=integer(dim(m)[1]))$clust
}

cluster.correlations <- function(matr,cutoff){
    d <- dim(matr);
    .C("clusterCorrelations",as.integer(d),as.numeric(t(matr)),as.numeric(cutoff),clust=integer(d[1]))$clust
}

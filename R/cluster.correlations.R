cluster.correlations <- function(matrix,cutoff){
    d <- dim(matrix);
    .C("clusterCorrelations",as.integer(d),as.numeric(matrix),as.numeric(cutoff),clust=integer(d[1]))$clust
}

cluster.correlations <- function(m,cutoff){
    d <- dim(m);
    .C("clusterCorrelations",as.integer(d),as.numeric(t(m)),as.numeric(cutoff),clust=integer(d[1]))$clust
}

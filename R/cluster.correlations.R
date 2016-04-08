cluster.correlations <- function(m,cutoff){
    d  <- dim(m);
    cl <- .C("clusterCorrelations",as.integer(d),as.numeric(t(m)),as.numeric(cutoff),clust=integer(d[1]))$clust
    df <- data.frame( cl=cl )
    df <- cbind(df,1:length(cl))
    colnames(df) <- c("cl","ind")
    clustering <- aggregate( df, by=list(Index=df$cl), function(x){x} )
    bigfirst  <- c(1)
    if( dim(clustering)[1] != 1 ){
        bigfirst <- order( unlist( lapply(clustering[,3],length) ), decreasing = T )
    }
    return(clustering[bigfirst,3])
}

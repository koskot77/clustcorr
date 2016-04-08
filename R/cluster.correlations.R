cluster.correlations <- function(m,cutoff){
    d  <- dim(m);
    cl <- .C("clusterCorrelations",as.integer(d),as.numeric(t(m)),as.numeric(cutoff),clust=integer(d[1]))$clust
    if( length( unique(cl) ) > 1 ){
        df <- data.frame( cl=cl )
        df <- cbind(df,1:length(cl))
        colnames(df) <- c("cl","ind")
        clustering   <- aggregate( df, by=list(Index=df$cl), function(x){x} )
        bigfirst     <- order( unlist( lapply(clustering[,3],length) ), decreasing = T )
        return(clustering[bigfirst,3])
    } else {
        # if all were merged in one cluster
        return(list('1'=list(cl)))
    }
}

recluster.correlations <- function(m,cutoff=NULL,nclust=NULL){
    if( !is.null(cutoff) ){
        if( !is.null(nclust) )  print("Both 'cutoff' and 'nclust' arguments are provided, ignoring 'nclust'")
        cl <- .C("reclusterCorrelationsAtCutoff",as.numeric(cutoff),clust=integer(dim(m)[1]))$clust
    } else {
        if( is.null(nclust) ){ print("Neither 'cutoff' nor 'nclust' arguments are provided, returning NULL"); return(NULL) }
        cl <- .C("reclusterCorrelationsForNclust",as.integer(nclust),clust=integer(dim(m)[1]))$clust
    }
    df <- data.frame( cl=cl )
    df <- cbind(df,1:length(cl))
    colnames(df) <- c("cl","ind")
    clustering <- aggregate( df, by=list(Index=df$cl), function(x){x} )
    bigfirst   <- order( unlist( lapply(clustering[,3],length) ), decreasing = T )
    return(clustering[bigfirst,3])
}

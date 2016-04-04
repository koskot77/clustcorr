cross.correlations <- function(m1,m2){
    d1 <- dim(m1)
    d2 <- dim(m2)
    if( d1[2] != d2[2] ){
        print("Windowing is not implemented")
        return(NULL)
    }
    gr <- .C("crossCorrelations",as.integer(d1),as.numeric(t(m1)),as.integer(d2),as.numeric(t(m2)),grouping=integer(d1[1]+d2[1]))$grouping
    gr1 <- data.frame( ind=seq(1:d1[1]), best.match=gr[1:d1[1]] )
    gr2 <- data.frame( int=seq(1:d2[1]), best.match=gr[d1[1]+1:d2[1]] )
    print("post-processing m1 -> m2")
    gr1 <- cbind(gr1, correlation = apply(gr1,1,function(x){ ifelse(x[2]!=0,cor(as.vector(m1[x[1],]),as.vector(m2[x[2],])),-10) }) )
    print("post-processing m2 -> m1")
    gr2 <- cbind(gr2, correlation = apply(gr2,1,function(x){ ifelse(x[2]!=0,cor(as.vector(m2[x[1],]),as.vector(m1[x[2],])),-10) }) )
    list('first2second'=gr1,'second2first'=gr2)
}

cross.correlations <- function(m1,m2,cutoff){
    d1 <- dim(m1);
    d2 <- dim(m2);
    gr <- .C("crossCorrelations",as.integer(d1),as.numeric(t(m1)),as.integer(d2),as.numeric(t(m2)),as.numeric(cutoff),grouping=integer(d1[1]+d2[1]))$grouping
    gr1 <- gr[1:d1[1]]
    gr2 <- gr[d1[1]+1:d2[1]]
    list('first2second'=gr1,'second2first'=gr2)
}

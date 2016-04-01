setcores <- function(n){
    .C("setCores",as.integer(n))
    return(NULL)
}

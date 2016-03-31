#include"utilities.h"
#include<math.h>

// mean, sd, and cor functions equivalent to their R analogs
double mean(const double *vals, size_t size){
    double sum = 0;
    for(unsigned int i=0; i<size; i++) sum += vals[i];
    return sum/size;
}

double sd(const double *vals, size_t size){
    double sum1 = 0, sum2 = 0;
    for(unsigned int i=0; i<size; i++){
        sum2 += vals[i]*vals[i];
        sum1 += vals[i];
    }
    return sqrt( (sum2-sum1*sum1/size)/(size-1) ); // unbiased
}

double cor(const double *a, double meanA, double sdA, const double *b, double meanB, double sdB, size_t size){
    double sum = 0;
    for(unsigned int i=0; i<size; i++)
        sum += a[i] * b[i];
    return size/double(size-1) * (sum/size-meanA*meanB) / sdA / sdB ;
}

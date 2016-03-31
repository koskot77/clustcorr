#ifndef utilities_h
#define utilities_h

// mean, sd, and cor functions equivalent to their R analogs
#include<stddef.h>

double mean(const double *vals, size_t size);
double sd  (const double *vals, size_t size);
double cor (const double *a, double meanA, double sdA, const double *b, double meanB, double sdB, size_t size);

#endif

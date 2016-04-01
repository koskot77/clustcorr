#include<string.h>
#include<stdio.h>
#include<math.h>
#include<thread>
#include<future>
#include"utilities.h"
#include"quicksort.h"
#include"unionfind.h"

// for stand-alone usage in R:
// g++ -Wl,--no-as-needed -g -Wall -std=c++11 -c -fPIC *.cc
// R CMD SHLIB -o correlator.so *.o -lstdc++ -lpthread
// R
// dyn.load("correlator.so")
// source("../R/cluster.correlations.R")

size_t nCores = 1;

// global copy of the pointer to input matrix
const double **source = 0;
size_t nRows = 0, nColumns = 0;

// global copy of the intermediate results (need to recompute when input changes)
double *cachedMean = 0, *cachedSd = 0;
float  *correlations = 0;        // use single precision to save some space
unsigned long long *indices = 0; // keep in mind that addressing scales N^2

extern "C" {

void cleanUp(void){
    // clean up if needed
    if( correlations ){ delete [] correlations; correlations = 0; }
    if( cachedMean   ){ delete [] cachedMean;   cachedMean   = 0; }
    if( cachedSd     ){ delete [] cachedSd;     cachedSd     = 0; }
    if( indices      ){ delete [] indices;      indices      = 0; }
}
void setCores(int *n){ nCores = *n; }
}

void computeBlockOfCorrelations(unsigned long long begin, unsigned long long end){
    const double (*series)[nColumns] = (const double (*)[nColumns])source;
    size_t length = nColumns;
    for(unsigned long long index=begin; index<end; index++){
        unsigned long row1  = (unsigned long)( (1+sqrt(1+8*index))/2 );
        unsigned long row2  = index - row1*(row1-1)/2;
        correlations[index] = -fabs( cor(series[row1],cachedMean[row1],cachedSd[row1],
                                         series[row2],cachedMean[row2],cachedSd[row2], length ) );
    }
}

void computeCorrelationsAndSort(void){
    cleanUp();

    // allocate space for all-with-all pair correlations 
    correlations        = new  float  [nRows*(nRows-1)/2];
    bzero(correlations, sizeof(float) *nRows*(nRows-1)/2);

    double (*series)[nColumns] = (double (*)[nColumns])(source);

    // caching means and sds
    cachedMean = new double [nRows];
    cachedSd   = new double [nRows];
    
    for(unsigned int row=0; row<nRows; row++){
        cachedMean[row] = mean(series[row],nColumns);
        cachedSd  [row] = sd  (series[row],nColumns);
    }

    // multicore calculations
    printf("Correlating all with all\n");
    const size_t maxNumThreads = nCores;
    std::future<void> results [ maxNumThreads ];
    const unsigned long long maxIndex  =  nRows*(nRows-1)/2;
    const unsigned int       maxBlocks = (nRows>100 ? 100 : 1);
    for(unsigned int block=0; block<maxBlocks; block++){
        unsigned long long begin =  block    * maxIndex / maxBlocks ;
        unsigned long long   end = (block+1) * maxIndex / maxBlocks ;

        // identify a free thread
        size_t freeThread = 0;
        for( ;  results[freeThread].valid() &&
                results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ; )
            if( freeThread == maxNumThreads-1 ) freeThread = 0; else freeThread++;

        // submit
        results[freeThread] = std::async(std::launch::async, computeBlockOfCorrelations, begin, end);
        printf("  finished  block %d/%d in thread %ld\n",block,maxBlocks,freeThread);
    }
    printf("  finalizing ... \n");

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

    /// Beginning of the sorting step

    // keep track of original indices while sorting
    indices = new unsigned long long [nRows*(nRows-1)/2];
    for(unsigned long long i=0; i<nRows*(nRows-1)/2; i++) indices[i] = i;

    printf("Sorting\n");
    unsigned long long size = nRows*(nRows-1)/2;
    splitQuickSort( correlations, indices, size );
}

extern "C" {

void clusterCorrelations(int *dim, double *inMatrix, double *cutoffDistance, int *outClustering){

    if( source  != (const double**)inMatrix || nRows != size_t(dim[0]) || nColumns != size_t(dim[1]) ){
        source   = (const double**)inMatrix;
        nRows    = dim[0];
        nColumns = dim[1];
        computeCorrelationsAndSort();
    }

    printf("Clustering\n");
    UnionFind uf(nRows);
    for(unsigned long long i=0; i<nRows*(nRows-1)/2; i++){
        unsigned long long index = indices[i];
        unsigned long row1  = (unsigned long)( (1+sqrt(1+8*index))/2 );
        unsigned long row2  = index - row1*(row1-1)/2;
        if( -correlations[i] < *cutoffDistance ) break;
        int node1 = row1+1;
        int node2 = row2+1;
        int cluster1 = uf.findCluster(node1);
        int cluster2 = uf.findCluster(node2);
        if( cluster1 != cluster2 ) uf.joinClusters(cluster1,cluster2);
    }

    for(unsigned int row=0; row<nRows; row++)
        outClustering[row] = uf.findCluster(row+1);

}

}

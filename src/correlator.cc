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

// global copy of the pointer to the input matrix
const double **source = 0;
size_t nRows = 0, nColumns = 0;
// global copy of the intermediate results (need to recompute when input changes)
double *cachedMean = 0, *cachedSd = 0;
float  *correlations = 0;        // use single precision to save some space
unsigned long long *indices = 0; // keep in mind that addressing scales N^2

// same for cross correlations between two different input matricies
const double **source1 = 0, **source2 = 0;
size_t nRows1 = 0, nColumns1 = 0;
size_t nRows2 = 0, nColumns2 = 0;
double *cachedMean1 = 0, *cachedSd1 = 0;
double *cachedMean2 = 0, *cachedSd2 = 0;

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

void reclusterCorrelations(double *cutoffDistance, int *outClustering){
    if( !source ){
        printf("Clustering not yet done, nothing to recluster\n");
        return;
    }

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

void computeBlockOfCrossCorrelations12(unsigned long long begin, unsigned long long end, int *grouping12){
    const double (*series1)[nColumns1] = (const double (*)[nColumns1])source1;
    const double (*series2)[nColumns2] = (const double (*)[nColumns2])source2;
    size_t length = nColumns1; // for now assumed to be = nCulumns2
    for(unsigned int row1=begin; row1<end; row1++){

        unsigned int bestMatchIndex = 0;
        double       bestMatchValue = 10;

        for(unsigned int row2=0; row2<nRows2; row2++){
            double correlation = -fabs( cor(series1[row1],cachedMean1[row1],cachedSd1[row1],
                                            series2[row2],cachedMean2[row2],cachedSd2[row2], length ) );
            if( correlation < bestMatchValue ){
                bestMatchValue = correlation;
                bestMatchIndex = row2;
            }
        }
        grouping12[row1] = bestMatchIndex + 1; // +1 to be consistent with R indexing
    }
}

void computeBlockOfCrossCorrelations21(unsigned long long begin, unsigned long long end, int *grouping21){
    const double (*series1)[nColumns1] = (const double (*)[nColumns1])source1;
    const double (*series2)[nColumns2] = (const double (*)[nColumns2])source2;
    size_t length = nColumns1; // for now assumed to be = nCulumns2
    for(unsigned int row2=begin; row2<end; row2++){

        unsigned int bestMatchIndex = 0;
        double       bestMatchValue = 10;

        for(unsigned int row1=0; row1<nRows1; row1++){
            double correlation = -fabs( cor(series1[row1],cachedMean1[row1],cachedSd1[row1],
                                            series2[row2],cachedMean2[row2],cachedSd2[row2], length ) );
            if( correlation < bestMatchValue ){
                bestMatchValue = correlation;
                bestMatchIndex = row1;
            }
        }

        grouping21[row2] = bestMatchIndex + 1; // +1 to be consistent with R indexing
    }
}


void computeBestCrossCorrelations(int *grouping12, int *grouping21){

    double (*series1)[nColumns1] = (double (*)[nColumns1])(source1);
    double (*series2)[nColumns2] = (double (*)[nColumns2])(source2);

    if( cachedMean1  ){ delete [] cachedMean1;  cachedMean1  = 0; }
    if( cachedSd1    ){ delete [] cachedSd1;    cachedSd1    = 0; }
    if( cachedMean2  ){ delete [] cachedMean2;  cachedMean2  = 0; }
    if( cachedSd2    ){ delete [] cachedSd2;    cachedSd2    = 0; }

    // caching means and sds
    cachedMean1 = new double [nRows1];
    cachedSd1   = new double [nRows1];
    cachedMean2 = new double [nRows2];
    cachedSd2   = new double [nRows2];
    
    for(unsigned int row=0; row<nRows1; row++){
        cachedMean1[row] = mean(series1[row],nColumns1);
        cachedSd1  [row] = sd  (series1[row],nColumns1);
    }

    for(unsigned int row=0; row<nRows2; row++){
        cachedMean2[row] = mean(series2[row],nColumns2);
        cachedSd2  [row] = sd  (series2[row],nColumns2);
    }

    // multicore calculations
    printf("Correlating all in one with all in another\n");
    const size_t maxNumThreads = nCores;
    std::future<void> results [ maxNumThreads ];

    // correlate first with the second
    unsigned long long maxIndex  =  nRows1;
    unsigned int       maxBlocks = (nRows1>100 ? 100 : 1);
    for(unsigned int block=0; block<maxBlocks; block++){
        unsigned long long begin =  block    * maxIndex / maxBlocks ;
        unsigned long long   end = (block+1) * maxIndex / maxBlocks ;

        // identify a free thread
        size_t freeThread = 0;
        for( ;  results[freeThread].valid() &&
                results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ; )
            if( freeThread == maxNumThreads-1 ) freeThread = 0; else freeThread++;

        // submit
        results[freeThread] = std::async(std::launch::async, computeBlockOfCrossCorrelations12, begin, end, grouping12);
        printf("  finished  block %d/%d in thread %ld\n",block,maxBlocks,freeThread);
    }
    printf("  finalizing12 ... \n");

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();


    // correlate second with the first 
    maxIndex  =  nRows2;
    maxBlocks = (nRows2>100 ? 100 : 1);
    for(unsigned int block=0; block<maxBlocks; block++){
        unsigned long long begin =  block    * maxIndex / maxBlocks ;
        unsigned long long   end = (block+1) * maxIndex / maxBlocks ;

        // identify a free thread
        size_t freeThread = 0;
        for( ;  results[freeThread].valid() &&
                results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ; )
            if( freeThread == maxNumThreads-1 ) freeThread = 0; else freeThread++;

        // submit
        results[freeThread] = std::async(std::launch::async, computeBlockOfCrossCorrelations21, begin, end, grouping21);
        printf("  finished  block %d/%d in thread %ld\n",block,maxBlocks,freeThread);
    }
    printf("  finalizing21 ... \n");

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

}

extern "C" {

void crossCorrelations(int *dim1, double *inMatrix1,int *dim2, double *inMatrix2, int *outGrouping){

    source1   = (const double**)inMatrix1;
    source2   = (const double**)inMatrix2;
    nRows1   = dim1[0];
    nColumns1= dim1[1];
    nRows2   = dim2[0];
    nColumns2= dim2[1];
    if( nColumns1 != nColumns2 ){
        printf("Windowing is not implemented\n");
        return;
    }

    int *outGrouping12 = outGrouping;
    int *outGrouping21 = outGrouping + nRows1;

    computeBestCrossCorrelations(outGrouping12,outGrouping21);

}

}

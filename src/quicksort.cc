#include "quicksort.h"
#include <stdio.h>
#include <thread>
#include <future>

extern size_t nCores;

// Canonical quick sort for ... sorting

enum {BEGIN_POS=1, MEDIAN_POS=2, END_POS=3};
unsigned int pivotMethod = MEDIAN_POS;

unsigned long long pivotElementPos(float *array, unsigned long long length){
    if( length<1 ) return 0;//exit(1);
//    return (rand()%length);
    switch( pivotMethod ){
        case BEGIN_POS  : return 0; break;
        case END_POS    : return length-1; break;
        case MEDIAN_POS :
            unsigned long long middlePos = ( length % 2 ? length/2 : length/2-1 );
            float a = array[0];
            float b = array[middlePos];
            float c = array[length-1];
            if( (a>=b && b>=c) || (a<=b && b<=c) ){ return middlePos; }
            if( (b>=a && a>=c) || (b<=a && a<=c) ){ return 0; }
            if( (a>=c && c>=b) || (a<=c && c<=b) ){ return length-1; }
            break;
    }
    return 0; //exit(1);
}

void swap(float &i,              float &j)             { float tmp=j; j=i; i=tmp; }
void swap(unsigned long long &i, unsigned long long &j){ int   tmp=j; j=i; i=tmp; }

unsigned long long partition(float *array, unsigned long long pivotPos, unsigned long long *payload, unsigned long long length){
    if( length <= 1 ) return 0;
    if( pivotPos >= length ) return 0; //exit(1);

    // preprocessing step
    float pivot = array[pivotPos];
    if( pivotPos != 0 ){
        swap( array  [pivotPos], array  [0] );
        swap( payload[pivotPos], payload[0] );
        pivotPos = 0;
    }

    unsigned long long i,j,k=1;
    for(i=1,j=1; j<length; j++){
        if( array[j] < pivot ){
            if( i != j ){
                swap( array  [i], array  [j] );
                swap( payload[i], payload[j] );
            }
            i++;
        }
        if( array[j] == pivot ) k++;
    }

    if( i == 1 && k == length ) // all elements are identical
        return i;

    swap( array  [i-1], array  [0] );
    swap( payload[i-1], payload[0] );

    return i;
}

void recure(float *array, unsigned long long *payload, unsigned long long length){

    unsigned long long i = partition(array, pivotElementPos(array,length), payload, length);

    float *first  = &(array[0]);
    float *second = &(array[i]);
    unsigned long long *plfirst  = &(payload[0]);
    unsigned long long *plsecond = &(payload[i]);

    if( i>1 )
        recure(first,  plfirst, i-1);
    if( length-i>1 )
        recure(second, plsecond, length-i);

}

void splitQuickSort(float *array, unsigned long long *payload, unsigned long long length){

    const unsigned int maxDepth = 5, maxBlocks = (1<<maxDepth);

    if( length < maxBlocks ){
        recure(array, payload, length);
        return;
    } 

    // consider those a perfectly balanced tree-like structure
    float              *subarray  [maxBlocks*2];
    unsigned long long *subpayload[maxBlocks*2];
    unsigned long long  sublength [maxBlocks*2];

    subarray  [0] = array;
    subpayload[0] = payload;
    sublength [0] = length;

    // unroll recursion into two linear loops
    for(unsigned int depth=0; depth<maxDepth; depth++){

        unsigned int subTasksSeen = (1<<depth)-1;

        for(unsigned int subtask=0; subtask<unsigned(1<<depth); subtask++){

            unsigned long long pos = partition( subarray  [ subtask + subTasksSeen ],
                                                pivotElementPos(
                                                    subarray [ subtask + subTasksSeen ],
                                                    sublength[ subtask + subTasksSeen ]
                                                ),
                                                subpayload[ subtask + subTasksSeen ],
                                                sublength [ subtask + subTasksSeen ] );

            unsigned int subTasksAtThisDepth = (1<<(depth+1))-1;

            subarray  [ 2*subtask+0 + subTasksAtThisDepth ] = subarray  [ subtask + subTasksSeen ];
            subarray  [ 2*subtask+1 + subTasksAtThisDepth ] = subarray  [ subtask + subTasksSeen ] + pos;
            subpayload[ 2*subtask+0 + subTasksAtThisDepth ] = subpayload[ subtask + subTasksSeen ];
            subpayload[ 2*subtask+1 + subTasksAtThisDepth ] = subpayload[ subtask + subTasksSeen ] + pos;
            sublength [ 2*subtask+0 + subTasksAtThisDepth ] = pos-1;
            sublength [ 2*subtask+1 + subTasksAtThisDepth ] = sublength [ subtask + subTasksSeen ] - pos;
        }
    }

    const size_t maxNumThreads = nCores;
    std::future<void> results [ maxNumThreads ];

    for(unsigned int block=0; block<maxBlocks; block++){

        unsigned int offset = maxBlocks - 1;

        if( sublength[block+offset] == 0 ) continue;

        // identify a free thread
        size_t freeThread = 0;
        for( ;  results[freeThread].valid() &&
                results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ; )
            if( freeThread == maxNumThreads-1 ) freeThread = 0; else freeThread++;

        // submit
        results[freeThread] = std::async(std::launch::async, recure, subarray[block+offset], subpayload[block+offset], sublength[block+offset]);

        printf("  finished  block %d/%d in thread %ld\n",block,maxBlocks,freeThread);
    }
    printf("  finalizing ... \n");

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

}

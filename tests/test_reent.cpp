
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <boost/thread.hpp>

#include "swmm5.h"

#define NUM_THREADS 2

using namespace std;

void swmmThread(long i)
{
    int errorcode = 0;
    SWMM_ProjectHandle ph;

    string prefix = "example_";
    string suffix = ".inp";
    string input = prefix + to_string(static_cast<long long>(i)) + suffix;

    suffix = ".rpt";
    string report = prefix + to_string(static_cast<long long>(i)) + suffix;

    suffix = ".out";
    string output = prefix + to_string(static_cast<long long>(i)) + suffix;

    printf("Thread #%ld starting SWMM ...\n", i);

    swmm_alloc_project(&ph);
    errorcode = swmm_run_project(ph, input.c_str(), report.c_str(), output.c_str());
    swmm_free_project(&ph);

    printf("Thread #%ld SWMM done. Status = %d\n", i, errorcode);
}

int main(int argc, char *argv[])
{
    long i;
    boost::thread *threads[NUM_THREADS];

    for (i = 0; i < NUM_THREADS; i++) {
        threads[i] = new boost::thread(swmmThread, i);
        printf("Main: creating thread %ld.\n", i);
    }

    for (i = 0; i < NUM_THREADS; i++) {
        threads[i]->join();
        printf("Main: joining thread %ld.\n", i);
        delete threads[i];
    }

    printf("Main: program completed. Exiting.\n");
    return(0);
}

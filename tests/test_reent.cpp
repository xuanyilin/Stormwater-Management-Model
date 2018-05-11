
#include <string>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "swmm5.h"

#define NUM_THREADS 2

using namespace std;

void *swmm_thread(void *threadid)
{
    long tid;
    SWMM_ProjectHandle ph;

    long errorcode = 0;
    tid = (intptr_t)threadid;

    string prefix = "example1_";
    string suffix = ".inp";
    string input = prefix + to_string(tid) + suffix;

    suffix = ".rpt";
    string report = prefix + to_string(tid) + suffix;

    suffix = ".out";
    string output = prefix + to_string(tid) + suffix;

    printf("Thread #%ld starting SWMM ...\n", tid);

    swmm_alloc_project(&ph);
    errorcode = swmm_run_project(ph, input.c_str(), report.c_str(), output.c_str());
    swmm_free_project(&ph);

    printf("Thread #%ld SWMM done. Status = %ld\n", tid, errorcode);
    pthread_exit((void *)(intptr_t)errorcode);

    return((void *)(intptr_t)errorcode);
}

int main(int argc, char *argv[])
{
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;

    int rc;
    long t;
    long status = 0;

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(t = 0; t < NUM_THREADS; t++){
        printf("In main: creating thread %ld\n", t);
        rc = pthread_create(&threads[t], &attr, swmm_thread, (void *)(intptr_t)t);
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);

    for(t = 0; t < NUM_THREADS; t++) {
        rc = pthread_join(threads[t], NULL); //(void **)(intptr_t)status);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }

        printf("Main: completed join with thread %ld having a status of %ld\n",
                t, status);
    }

    printf("Main: program completed. Exiting.\n");
    pthread_exit(NULL);
    return(0);
}

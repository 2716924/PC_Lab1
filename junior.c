// keep in mind that this is false sharing


#include <stdio.h>
#include <omp.h>
static long num_steps = 1000000;
double step;
#define num_threads 2
int main(){
    int n_threads;
    double pi, sum[num_threads]; //array of doubles to store sum of each thread
    step = 1.0/(double)num_steps;
    omp_set_num_threads(num_threads); //does not really make sense but hitaku yini
    #pragma omp parallel 
    {
        int i,id,tthreads;
        double x = 0.0;
        id = omp_get_thread_num();
        tthreads = omp_get_num_threads();
        if (id==0) n_threads = tthreads;
        //we now populate the sum array
        for(i=id,sum[id]=0.0;i<num_steps;i+=tthreads){
            x = (i+0.5)*step;
            sum[id] = sum[id] + 4.0/(1.0+x*x);
        }

    }
    for(int i=0;i<n_threads;i++){
        pi+= step*sum[i];
    }
    printf("pi" = %f\n",pi); 
    return 0;
}


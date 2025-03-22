#include <stdio.h>
#include <omp.h>
#include <stdbool.h>

static long num_steps = 1000000;

int main(){
    int i,tthreads,id;
    double step,pi,sum=0.0,x;
    step = 1.0/(double)num_steps;


        #pragma omp parallel 
        {
            tthreads = omp_get_num_threads();
            id = omp_get_thread_num();
            for(i=id,sum=0.0;i<num_steps;i+=tthreads){
                x = (i+0.5)*step;
                sum = sum + 4.0/(1.0+x*x);
            }
        }


        pi = step*sum;
        printf("pi = %f\n",pi);
       
    
    

    return 0;

}

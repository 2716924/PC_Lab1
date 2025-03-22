#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double compute_pi(double step);
double false_sharing(double step);
double method_2(double step);
double parallel_pi(double step);

void usage(char prog_name[]);

static long num_steps = 1000000;
#define num_threads 12  // Set the number of threads

int main(int argc, char **argv) {
    double step = 1.0 / (double)num_steps;
    double start_time, run_time;
    double pi_seq, pi_false_sharing, pi_race_condition, pi_parallel;
    double time_seq, time_false_sharing, time_race_condition, time_parallel;

    // Sequential computation
    start_time = omp_get_wtime();
    pi_seq = compute_pi(step);
    time_seq = omp_get_wtime() - start_time;

    // Parallel with false sharing
    start_time = omp_get_wtime();
    pi_false_sharing = false_sharing(step);
    time_false_sharing = omp_get_wtime() - start_time;

    // Parallel with race condition
    start_time = omp_get_wtime();
    pi_race_condition = method_2(step);
    time_race_condition = omp_get_wtime() - start_time;

    // Parallel with no race condition
    start_time = omp_get_wtime();
    pi_parallel = parallel_pi(step);
    time_parallel = omp_get_wtime() - start_time;

    // Print results
    printf("Sequential: pi with %ld steps is %.6f in %.6f seconds\n", num_steps, pi_seq, time_seq);
    printf("Parallel with false sharing: pi with %ld steps is %.6f in %.6f seconds using %d threads\n",
           num_steps, pi_false_sharing, time_false_sharing, num_threads);
    printf("Speedup: %.6f\n", time_seq / time_false_sharing);

    printf("Parallel with race condition: pi with %ld steps is %.6f in %.6f seconds using %d threads\n",
           num_steps, pi_race_condition, time_race_condition, num_threads);
    printf("Speedup: %.6f\n", time_seq / time_race_condition);

    printf("Parallel with no race and parallel: pi with %ld steps is %.6f in %.6f seconds using %d threads\n",
           num_steps, pi_parallel, time_parallel, num_threads);
    printf("Speedup: %.6f\n", time_seq / time_parallel);

    return EXIT_SUCCESS;
}

/*--------------------------------------------------------------------
 * Function:    compute_pi
 * Purpose:     Compute number pi in serial
 * In arg:      step
 */  
double compute_pi(double step) {
	double pi, x, sum = 0.0;
	for (int i = 0; i < num_steps; i++) {
		x = (i + 0.5) * step;
		sum += 4.0 / (1.0 + x * x);
	}
	pi = step * sum;
	return pi;
}

/*--------------------------------------------------------------------
 * Function:    usage
 * Purpose:     Print command line for function
 * In arg:      prog_name
 */
void usage(char prog_name[]) {
   fprintf(stderr, "usage:  %s <number of times to run>\n", prog_name);
}

/*--------------------------------------------------------------------
 * Function:    method_2
 * Purpose:     Compute pi with intentional false sharing using OpenMP
 * In arg:      step
 */
double method_2(double step){
	double pi = 0.0;
	double sum = 0.0;
	double x;
	#pragma omp parallel 
	{
		int id = omp_get_thread_num();
		int tthreads = omp_get_num_threads();
		for(int i=id; i<num_steps; i+=tthreads){
			x = (i + 0.5) * step;
			sum += 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;
	return pi;
}

/*--------------------------------------------------------------------
 * Function:    false_sharing
 * Purpose:     Compute pi with intentional false sharing using OpenMP
 * In arg:      step
 */
double false_sharing(double step) {
    double pi = 0.0;
    double sum = 0.0;

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        double x, local_sum = 0.0; // Use local sum

        for (int i = id; i < num_steps; i += nthrds) {
            x = (i + 0.5) * step;
            local_sum += 4.0 / (1.0 + x * x);
        }

        
        sum += local_sum;  // Safely update sum
    }

    pi = step * sum;
    return pi;
}

/*--------------------------------------------------------------------
 * Function:    parallel_pi
 * Purpose:     Compute pi in a parallel and thread-safe manner
 * In arg:      step
 */
double parallel_pi(double step) {
    double pi = 0.0;
    
    omp_set_num_threads(num_threads);  // Use correct `num_threads`

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        double x, local_sum = 0.0;

        for (int i = id; i < num_steps; i += nthrds) {
            x = (i + 0.5) * step;
            local_sum += 4.0 / (1.0 + x * x);
        }

        #pragma omp atomic
        pi += local_sum * step;  // Safely update global `pi`
    }

    return pi;
}

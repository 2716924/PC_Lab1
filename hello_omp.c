#include <stdio.h>
#include <omp.h>
int main()
{
	//Parallel region with default number of threads
	//int num;
	//num = omp_get_num_procs();
	//omp_set_num_threads(4);
	#pragma omp parallel //num_threads(num)
	//Start of the parallel region
	{
		//Runtime library function to return a thread number
		int ID = omp_get_thread_num();
		printf("Hello World! (Thread %d)\n", ID);
	}//End of the parallel region
    return 0;
}

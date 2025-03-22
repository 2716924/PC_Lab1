echo ---------hello_omp is started----------
./hello_omp |& tee -a terminal_output.out
echo ---------hello_omp is done-------------
echo
echo ---------parallel_region is started-----
./parallel_region |& tee -a terminal_output.out
echo ---------parallel_region is done--------
echo
echo -----------pi is started----------------
./pi 3 |& tee -a terminal_output.out
echo -----------pi is done--------------------

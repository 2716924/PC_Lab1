#This is a comment
P1=pi
P2=parallel_region
P3=hello_omp
#INC="./inc"
#FLAGS=-I$(INC)
OMPFLAG=-fopenmp
CC=gcc
CFLAGS = -g -Wall

all: hello_omp parallel_region $(P1) $(P2) $(P3)

$(P1): $(P1).c
	$(CC) $(CFLAGS) $(OMPFLAG) $(P1).c -o $(P1)
$(P2): $(P2).c
	$(CC) $(CFLAGS) $(OMPFLAG) $(P2).c -o $(P2)
$(P3): $(P3).c
	$(CC) $(CFLAGS) $(OMPFLAG) $(P3).c -o $(P3)


clean:
	rm -vf $(P1) $(P2) $(P3)

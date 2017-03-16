CFLAGS = -lm -gstabs+

main: main_mpi.c 
	mpicc -ggdb -Wall -o main main_mpi.c Lab4_IO.c $(CFLAGS)

main_serial: main_serial.c 
	gcc -ggdb -Wall -o main_serial main_serial.c Lab4_IO.c $(CFLAGS)

data: datatrim.c
	gcc datatrim.c Lab4_IO.c -o datatrim $(CFLAGS)

serialtester: serialtester.c
	gcc serialtester.c Lab4_IO.c -o serialtester -lm

	



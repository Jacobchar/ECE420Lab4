CFLAGS = -lm

main: main.c 
	gcc -Wall -o main main.c Lab4_IO.c $(CFLAGS)

data: datatrim.c
	gcc datatrim.c Lab4_IO.c -o data_input $(CFLAGS)

serialtester: serialtester.c
	gcc serialtester.c Lab4_IO.c -o serialtester -lm

	



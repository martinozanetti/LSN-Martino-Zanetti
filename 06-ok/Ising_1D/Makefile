CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : ../../random/random.cpp ../../random/random.h
	$(CC) -c ../../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe
run : Monte_Carlo_ISING_1D.exe
	./Monte_Carlo_ISING_1D.exe

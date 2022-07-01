CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE.exe : MD.o random.o
	$(CC) random.o MD.o -o NVE.exe
MD.o : MD.cpp
	$(CC) -c MD.cpp -o MD.o $(CFLAGS)
random.o : ../../random/random.cpp ../../random/random.h
	$(CC) -c ../../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe

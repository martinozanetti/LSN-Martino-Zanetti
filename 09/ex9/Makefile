CC = g++
CFLAGS = -Wall -O3 --std=c++11


exe : main.o random.o mylib.o tsplib.o
	$(CC) random.o mylib.o tsplib.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : ../../random/random.cpp ../../random/random.h
	$(CC) -c ../../random/random.cpp -o random.o $(CFLAGS)
mylib.o : ../../mylib/mylib.cpp ../../mylib/mylib.h 
	$(CC) -c ../../mylib/mylib.cpp -o mylib.o $(CFLAGS)
tsplib.o : ../../tsplib/tsplib.cpp ../../tsplib/tsplib.h ../../tsplib/parameters.h
	$(CC) -c ../../tsplib/tsplib.cpp -o tsplib.o $(CFLAGS)
clean :
	rm *.o main.exe
clean-out :
	rm imgs/*.png cit/*.out bestLen/*.out
clean-all : clean clean-out
run : exe
	./main.exe
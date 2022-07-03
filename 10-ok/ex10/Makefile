## Nota: per qualche motivo anaconda si è messa in mezzo nell'installazione di mpich 
#  e scrivendo solo mpicxx il file veniva cercato nel bin di anaconda3. Ho dovuto reindirizzare
#  al percorso giusto /bin/mpicxx.
CXX = /bin/mpicxx
CFLAGS = -Wall -O3 --std=c++11


exe : main.o random.o mylib.o tsplib.o
	$(CXX) random.o mylib.o tsplib.o main.o -o main.exe
main.o : main.cpp
	$(CXX) -c main.cpp -o main.o $(CFLAGS)
random.o : ../../random/random.cpp ../../random/random.h
	$(CXX) -c ../../random/random.cpp -o random.o $(CFLAGS)
mylib.o : ../../mylib/mylib.cpp ../../mylib/mylib.h 
	$(CXX) -c ../../mylib/mylib.cpp -o mylib.o $(CFLAGS)
tsplib.o : ../../tsplib/tsplib.cpp ../../tsplib/tsplib.h ../../tsplib/parameters.h
	$(CXX) -c ../../tsplib/tsplib.cpp -o tsplib.o $(CFLAGS)
clean :
	rm *.o main.exe
clean-out :
	rm imgs/*.png cit/*.out bestLen/*.out
clean-all : clean clean-out
run : exe
	mpiexec -np 3 ./main.exe
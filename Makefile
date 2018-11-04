CC:=/usr/local/bin/g++-8 #update to your local compiler with OPENMP
OMP:=-fopenmp #update according to compiler

main : main.o fsm.o coarseinitial.o coarseupdate.o coarseupdatenewtheta.o fineupdate.o causalsweeping.o
	$(CC) $(OMP) -O3 -o main main.o coarseinitial.o fsm.o coarseupdate.o coarseupdatenewtheta.o fineupdate.o causalsweeping.o

main.o : main.cpp sweepdefs.h
	$(CC) $(OMP) -c main.cpp
fsm.o : fsm.cpp sweepdefs.h
	$(CC) $(OMP) -c fsm.cpp
coarseinitial.o : coarseinitial.cpp sweepdefs.h
	$(CC) $(OMP) -c coarseinitial.cpp
coarseupdate.o : coarseupdate.cpp sweepdefs.h
	$(CC) $(OMP) -c coarseupdate.cpp
coarseupdatenewtheta.o : coarseupdatenewtheta.cpp sweepdefs.h
	$(CC) $(OMP) -c coarseupdatenewtheta.cpp
fineupdate.o : fineupdate.cpp sweepdefs.h
	$(CC) $(OMP) -c fineupdate.cpp
causalsweeping.o : causalsweeping.cpp sweepdefs.h
	$(CC) $(OMP) -c causalsweeping.cpp

clean: 
	rm main.o coarseinitial.o fsm.o coarseupdate.o coarseupdatenewtheta.o fineupdate.o causalsweeping.o

EXEC=Indenteur
CC=g++

S=GRAPH

BITS=$(getconf LONG_BIT)
ifeq ($(BITS),64)
	FBITS=-m64
else
	FBITS=
endif

FLAGS=$(FBITS) -O3

all : Factor_Graph

Factor_Graph : main.o problems.o graph.o 
	$(CC) $(FLAGS) -o bin/Factor_Graph main.o problems.o graph.o 

main.o : $(S)/main.cpp $(S)/graph.h $(S)/problems.h
	$(CC) -c $(FLAGS) $(S)/main.cpp

problems.o : $(S)/problems.cpp $(S)/problems.h $(S)/graph.h
	$(CC) -c $(FLAGS) $(S)/problems.cpp

graph.o : $(S)/graph.cpp $(S)/graph.h 
	$(CC) -c $(FLAGS) $(S)/graph.cpp

#per cancellare i file oggetto fai make clean
clean :
	@rm *.o

rmexec :
	@rm $(EXEC)
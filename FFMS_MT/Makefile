#TARGET = cplex greedy local_search metaheuristic hybrid_metaheuristic
TARGET = greedy greedyClasses pruebas p-greedy GENETIC DERIVATIVE BARRAKUDA p-cmsa
CXXFLAGS = -ansi -O3 -fpermissive -std=c++17 

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /home/eileen/Desktop/CPLEX_Studio201/cplex
CONCERTDIR    = /home/eileen/Desktop/CPLEX_Studio201/concert
GCC = gcc
CCC = g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -fpermissive -w
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -lpthread -ldl
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread -lpthread -ldl
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
OBJS = Timer.o
all: ${TARGET}

#cplex: cplex.o
#	$(CCC) $(CCFLAGS) cplex.o -o cplex $(CCLNFLAGS)

#cplex.o: cplex.cpp
#	$(CCC) -c $(CCFLAGS) cplex.cpp -o cplex.o 

#hybrid_metaheuristic: hybrid_metaheuristic.o 
#	$(CCC) $(CCFLAGS) hybrid_metaheuristic.o -o hybrid_metaheuristic $(CCLNFLAGS)

#hybrid_metaheuristic.o: hybrid_metaheuristic.cpp
#	$(CCC) -c $(CCFLAGS) hybrid_metaheuristic.cpp -o hybrid_metaheuristic.o 
 
greedy: greedy.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

p-greedy: p-greedy.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

greedyClasses: greedyClasses.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

pruebas: pruebas.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

#derivativeAnnealing: derivativeAnnealing.cpp
#	${CCC} ${CXXFLAGS} -o $@ $^

GENETIC: GENETIC.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

DERIVATIVE: derivativeAnnealing.cpp
	${CCC} ${CXXFLAGS} -o $@ $^
#local_search: local_search.cpp
#	${CCC} ${CXXFLAGS} -o $@ $^

#metaheuristic: metaheuristic.cpp
#	${CCC} ${CXXFLAGS} -o $@ $^

p-cmsa: p-cmsa.o $(OBJS)
	$(CCC) $(CCFLAGS) p-cmsa.o $(OBJS) -o p-cmsa $(CCLNFLAGS)

p-cmsa.o: p-cmsa.cpp 
	$(CCC) -c $(CCFLAGS) p-cmsa.cpp -o p-cmsa.o 

BARRAKUDA: BARRAKUDA.o $(OBJS)
	$(CCC) $(CCFLAGS) BARRAKUDA.o $(OBJS) -o BARRAKUDA $(CCLNFLAGS)

BARRAKUDA.o: BARRAKUDA.cpp 
	$(CCC) -c $(CCFLAGS) BARRAKUDA.cpp -o BARRAKUDA.o 

clean:
	@rm -f *~ *.o ${TARGET} core


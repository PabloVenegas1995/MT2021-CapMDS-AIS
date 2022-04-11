TARGET = cplex greedy local_search metaheuristic hybrid_metaheuristic
CXXFLAGS = -ansi -O3 -fpermissive -std=c++17 

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /home/eileen/Desktop/Cplex_FFMS/CPLEX_Studio201/cplex/
CONCERTDIR    = /home/eileen/Desktop/Cplex_FFMS/CPLEX_Studio201/concert
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

all: ${TARGET}

cplex: cplex.o
	$(CCC) $(CCFLAGS) cplex.o -o cplex $(CCLNFLAGS)

cplex.o: cplex.cpp
	$(CCC) -c $(CCFLAGS) cplex.cpp -o cplex.o 

hybrid_metaheuristic: hybrid_metaheuristic.o 
	$(CCC) $(CCFLAGS) hybrid_metaheuristic.o -o hybrid_metaheuristic $(CCLNFLAGS)

hybrid_metaheuristic.o: hybrid_metaheuristic.cpp
	$(CCC) -c $(CCFLAGS) hybrid_metaheuristic.cpp -o hybrid_metaheuristic.o 
 
greedy: greedy.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

local_search: local_search.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

metaheuristic: metaheuristic.cpp
	${CCC} ${CXXFLAGS} -o $@ $^

clean:
	@rm -f *~ *.o ${TARGET} core


TARGET = barrakuda_solver_simple_cleanModel
#TARGET = CMSA_capmds_PRE-BCN-2019-t CMSA_capmds_PRE-BCN-2019 CMSA_capmds_Barcelona ilp_capmds_Barcelona CMSA_capmds_Barcelona_ModeloPersistente CMSA_capmds_Barcelona_ModeloPersistente-t CMSA_capmds_BCN_MP_Dinamico ilp_capmds_Barcelona_Barrakuda
#TARGET = ilp_capmds ilp_capmds_eficiente ilp_capmds_var CMSA_capmds CMSA_capmds_tuning CMSA_capmds_dep CMSA_capmds_dep_Blum CMSA_capmds_Barcelona ilp_capmds_Barcelona CMSA_capmds_Barcelona_ModeloPersistente CMSA_capmds_Barcelona_ModeloPersistente-t CMSA_capmds_BCN_MP_Dinamico
CXXFLAGS = -ansi -O3 -std=c++11 
#CXXFLAGS = -std=c++11 
#CXXFLAGS = -std=c++17 
CXXFLAGS-B = -ansi -O3 -fpermissive -std=c++17 -lstdc++fs
#CXXFLAGS = -std=gnu++11 
#CXXFLAGS = -std=c++0x
OBJS = Timer.o Random.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /home/eileen/Desktop/GitHubRepositories/Barrakuda/CPLEX_Studio128/cplex
CONCERTDIR    = /home/eileen/Desktop/GitHubRepositories/Barrakuda/CPLEX_Studio128/concert
#CCC = g++-5
CCC = g++

#CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++11 -fpermissive -w
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -fpermissive -w
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -lpthread -ldl -lstdc++fs
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread -lpthread -ldl
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
all: ${TARGET}

ilp_capmds: ilp_capmds.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds.o $(OBJS) -o ilp_capmds $(CCLNFLAGS)

ilp_capmds.o: ilp_capmds.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds.cpp -o ilp_capmds.o 


ilp_capmds_Barcelona: ilp_capmds_Barcelona.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds_Barcelona.o $(OBJS) -o ilp_capmds_Barcelona $(CCLNFLAGS)

ilp_capmds_Barcelona.o: ilp_capmds_Barcelona.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds_Barcelona.cpp -o ilp_capmds_Barcelona.o 

ilp_capmds_Barcelona_Barrakuda: ilp_capmds_Barcelona_Barrakuda.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds_Barcelona_Barrakuda.o $(OBJS) -o ilp_capmds_Barcelona_Barrakuda $(CCLNFLAGS)

ilp_capmds_Barcelona_Barrakuda.o: ilp_capmds_Barcelona_Barrakuda.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds_Barcelona_Barrakuda.cpp -o ilp_capmds_Barcelona_Barrakuda.o 

ilp_capmds_Final: ilp_capmds_Final.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds_Final.o $(OBJS) -o ilp_capmds_Final $(CCLNFLAGS)

ilp_capmds_Final.o: ilp_capmds_Final.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds_Final.cpp -o ilp_capmds_Final.o 

ilp_capmds_var: ilp_capmds_var.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds_var.o $(OBJS) -o ilp_capmds_var $(CCLNFLAGS)

ilp_capmds_var.o: ilp_capmds_var.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds_var.cpp -o ilp_capmds_var.o 

ilp_capmds_eficiente: ilp_capmds_eficiente.o $(OBJS)
	$(CCC) $(CCFLAGS) ilp_capmds_eficiente.o $(OBJS) -o ilp_capmds_eficiente $(CCLNFLAGS)

ilp_capmds_eficiente.o: ilp_capmds_eficiente.cpp
	$(CCC) -c $(CCFLAGS) ilp_capmds_eficiente.cpp -o ilp_capmds_eficiente.o 

CMSA_capmds: CMSA_capmds.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds.o $(OBJS) -o CMSA_capmds $(CCLNFLAGS)

CMSA_capmds.o: CMSA_capmds.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds.cpp -o CMSA_capmds.o 


CMSA_capmds_tuning: CMSA_capmds_tuning.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_tuning.o $(OBJS) -o CMSA_capmds_tuning $(CCLNFLAGS)

CMSA_capmds_tuning.o: CMSA_capmds_tuning.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_tuning.cpp -o CMSA_capmds_tuning.o 


CMSA_capmds_dep: CMSA_capmds_dep.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_dep.o $(OBJS) -o CMSA_capmds_dep $(CCLNFLAGS)

CMSA_capmds_dep.o: CMSA_capmds_dep.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_dep.cpp -o CMSA_capmds_dep.o 


CMSA_capmds_dep_Blum: CMSA_capmds_dep_Blum.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_dep_Blum.o $(OBJS) -o CMSA_capmds_dep_Blum $(CCLNFLAGS)

CMSA_capmds_dep_Blum.o: CMSA_capmds_dep_Blum.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_dep_Blum.cpp -o CMSA_capmds_dep_Blum.o 


CMSA_capmds_Barcelona: CMSA_capmds_Barcelona.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_Barcelona.o $(OBJS) -o CMSA_capmds_Barcelona $(CCLNFLAGS)

CMSA_capmds_Barcelona.o: CMSA_capmds_Barcelona.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_Barcelona.cpp -o CMSA_capmds_Barcelona.o 



CMSA_capmds_Barcelona_ModeloPersistente-t: CMSA_capmds_Barcelona_ModeloPersistente-t.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_Barcelona_ModeloPersistente-t.o $(OBJS) -o CMSA_capmds_Barcelona_ModeloPersistente-t $(CCLNFLAGS)

CMSA_capmds_Barcelona_ModeloPersistente-t.o: CMSA_capmds_Barcelona_ModeloPersistente-t.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_Barcelona_ModeloPersistente-t.cpp -o CMSA_capmds_Barcelona_ModeloPersistente-t.o 


CMSA_capmds_Barcelona_ModeloPersistente: CMSA_capmds_Barcelona_ModeloPersistente.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_Barcelona_ModeloPersistente.o $(OBJS) -o CMSA_capmds_Barcelona_ModeloPersistente $(CCLNFLAGS)

CMSA_capmds_Barcelona_ModeloPersistente.o: CMSA_capmds_Barcelona_ModeloPersistente.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_Barcelona_ModeloPersistente.cpp -o CMSA_capmds_Barcelona_ModeloPersistente.o 

CMSA_capmds_BCN_MP_Dinamico: CMSA_capmds_BCN_MP_Dinamico.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_BCN_MP_Dinamico.o $(OBJS) -o CMSA_capmds_BCN_MP_Dinamico $(CCLNFLAGS)

CMSA_capmds_BCN_MP_Dinamico.o: CMSA_capmds_BCN_MP_Dinamico.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_BCN_MP_Dinamico.cpp -o CMSA_capmds_BCN_MP_Dinamico.o 


CMSA_capmds_PRE-BCN-2019: CMSA_capmds_PRE-BCN-2019.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_PRE-BCN-2019.o $(OBJS) -o CMSA_capmds_PRE-BCN-2019 $(CCLNFLAGS)

CMSA_capmds_PRE-BCN-2019.o: CMSA_capmds_PRE-BCN-2019.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_PRE-BCN-2019.cpp -o CMSA_capmds_PRE-BCN-2019.o 


CMSA_capmds_PRE-BCN-2019-t: CMSA_capmds_PRE-BCN-2019-t.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA_capmds_PRE-BCN-2019-t.o $(OBJS) -o CMSA_capmds_PRE-BCN-2019-t $(CCLNFLAGS)

CMSA_capmds_PRE-BCN-2019-t.o: CMSA_capmds_PRE-BCN-2019-t.cpp
	$(CCC) -c $(CCFLAGS) CMSA_capmds_PRE-BCN-2019-t.cpp -o CMSA_capmds_PRE-BCN-2019-t.o 


brkga: brkga.cpp $(OBJS)
	${CXX} ${CXXFLAGS-B} $^ -o $@ -lstdc++fs

brkga_v2: brkga_v2.cpp $(OBJS)
	${CXX} ${CXXFLAGS-B} $^ -o $@ 

barrakuda_solver: barrakuda_solver.o $(OBJS)
	$(CXX) $(CXXFLAGS-B) $(CCFLAGS) barrakuda_solver.o $(OBJS) -o barrakuda_solver $(CCLNFLAGS) 

barrakuda_solver.o: barrakuda_solver.cpp
	$(CXX) -c $(CXXFLAGS-B)  $(CCFLAGS) barrakuda_solver.cpp  -o  barrakuda_solver.o 


barrakuda_solver_simple: barrakuda_solver_simple.o $(OBJS)
	$(CXX) $(CXXFLAGS-B) $(CCFLAGS) barrakuda_solver_simple.o $(OBJS) -o barrakuda_solver_simple $(CCLNFLAGS) 

barrakuda_solver_simple.o: barrakuda_solver_simple.cpp
	$(CXX) -c $(CXXFLAGS-B)  $(CCFLAGS) barrakuda_solver_simple.cpp  -o  barrakuda_solver_simple.o 


barrakuda_solver_simple_cleanModel: barrakuda_solver_simple_cleanModel.o $(OBJS)
	$(CXX) $(CXXFLAGS-B) $(CCFLAGS) barrakuda_solver_simple_cleanModel.o $(OBJS) -o barrakuda_solver_simple_cleanModel $(CCLNFLAGS) 

barrakuda_solver_simple_cleanModel.o: barrakuda_solver_simple_cleanModel.cpp
	$(CXX) -c $(CXXFLAGS-B)  $(CCFLAGS) barrakuda_solver_simple_cleanModel.cpp  -o  barrakuda_solver_simple_cleanModel.o 

clean:
	@rm -f *~ *.o ${TARGET} core


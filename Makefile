CXX = g++ 

N_STATE = 2
N_NOISE = 2
T_ORDER = 2

CURRENT_DIR = $(shell pwd)
EXTRA_ARGS = -DUSE_MAF1 -DVERBOSE -DN_STATE=$(N_STATE) -DN_NOISE=$(N_NOISE) -DT_ORDER=$(T_ORDER)

CXXFLAGS = -std=c++11 $(EXTRA_ARGS) -Wall -I. -I$(CURRENT_DIR)/aa

SOURCES_utils = ode_integr.cpp ode_dyn.cpp aa/aa.cpp aa/aa_mod2.cpp aa/aa_mod.cpp aa/interval.cpp
Include_utils = config_outer.h ode_integr.h ode_dyn.h aa/aa_mod.h aa/aa.h aa/aa_mod2.h aa/config.h aa/interval.h

SOURCES = $(SOURCES_utils) main.cpp

OBJECTS_utils = $(SOURCES_utils:%.cpp=%.o)

OBJECTS = $(SOURCES:%.cpp=%.o)

all: $(OBJECTS) main

main : $(SOURCES) $(OBJECTS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS_utils) $@.o

clean:
	-rm *.o aa/*.o main *.out

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

%.o: %.cpp $(Include_utils)
	$(CXX) $(CXXFLAGS) -c $(@:%.o=%.cpp) -o $@ 
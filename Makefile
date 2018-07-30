CXX = g++ 

N_STATE = 2
N_NOISE = 2
T_ORDER = 2

CURRENT_DIR = $(shell pwd)

CXXFLAGS = -Wall -DN_STATE=$(N_STATE) -DN_NOISE=$(N_NOISE) -DT_ORDER=$(T_ORDER) -I. -I$(CURRENT_DIR)/aa

SOURCES_utils = ode_integr.cpp ode_dyn.cpp aa/aa.cpp aa/interval.cpp

SOURCES = $(SOURCES_utils) main.cpp

OBJECTS_utils = $(SOURCES_utils:%.cpp=%.o)

OBJECTS = $(SOURCES:%.cpp=%.o)

all: $(OBJECTS) main

main : $(SOURCES) $(OBJECTS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS_utils) $@.o

clean:
	-rm *.o main

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c $(@:%.o=%.cpp) -o $@ 
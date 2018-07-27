LIB=libaf.a

AR = ar
ARFLAGS = -r

CC=g++

CFLAGS=-c -Wall
# EXTRA_ARGS= -DFAST_MULT -DFAST_DIV -DN_NOISE=3
EXTRA_ARGS= -DN_NOISE=3

LDFLAGS=

RM = rm -f 

SOURCES=aa.cpp interval.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=test
OBJECTS_BIN= $(OBJECTS) $(EXECUTABLE).o

static:	$(LIB)

$(LIB): $(SOURCES) $(OBJECTS)
		@echo "\n - Compiling $(LIB) -\n"
		@rm -f $(LIB)
		$(AR) $(ARFLAGS) $(EXTRA_ARGS) $(LIB) $(OBJECTS)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS_BIN) 
	$(CC) $(LDFLAGS) $(EXTRA_ARGS) $(OBJECTS_BIN) -o $@

.cpp.o:
	$(CC) $(EXTRA_ARGS) $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJECTS_BIN) $(LIB) $(EXECUTABLE)
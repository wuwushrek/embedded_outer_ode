LIB=libaf.a

AR = ar
ARFLAGS = -r

CC=g++

CFLAGS=-c -Wall

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
		$(AR) $(ARFLAGS) $(LIB) $(OBJECTS)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS_BIN) 
	$(CC) $(LDFLAGS) $(OBJECTS_BIN) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJECTS_BIN) $(LIB) $(EXECUTABLE)
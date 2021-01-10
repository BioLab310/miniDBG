EXECUTABLE=miniDBG

CXX = g++

CFLAGS = -lpthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
endif
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/opt/local/include
endif

$(EXECUTABLE): $(objects)
	$(CXX) -o $@ $(objects) $(CFLAGS)

objects = Main.o Basic.o Binary_Search.o Builder.o Counters.o Locate.o Read.o
	
all:$(objects)

$(objects): %.o :%.cpp
	$(CXX) -c $< -o $@

.PHONY : clean

clean:
	rm *.o $(EXECUTABLE)

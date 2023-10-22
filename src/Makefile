EXECUTABLE=MiniDBG

CXX = g++

CFLAGS = -lpthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
endif
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/opt/local/include
endif


objects = Main.o Basic.o Binary_Search.o Builder.o Counters.o Locate.o Read.o
	
all:$(objects) $(EXECUTABLE)

$(objects): %.o :%.cpp
	$(CXX) -c $< -o $@

$(EXECUTABLE): $(objects)
	$(CXX) -o $@ $(objects) $(CFLAGS)

.PHONY : clean

clean:
	@if [ -f "$(EXECUTABLE)" ]; then rm -f $(EXECUTABLE); fi
	@if [ -n "$(objects)" ]; then rm -f $(objects); fi

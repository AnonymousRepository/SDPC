CXX=g++
RM=rm -f
CXXFLAGS=-std=c++11

SRCS=main.cpp DPC.cpp
OBJS=main.o DPC.o

all: DPC_Express clean

DPC_Express: $(OBJS)
	$(CXX) $(CXXFLAGS) -o DPC_Express $(OBJS)

DPC.o: DPC.cpp DPC.h
	$(CXX) $(CXXFLAGS) -c DPC.cpp DPC.h
main.o: main.cpp DPC.h
	$(CXX) $(CXXFLAGS) -c main.cpp DPC.h
clean:
	$(RM) $(OBJS) *.h.gch
distclean: clean
	$(RM) DPC_Express

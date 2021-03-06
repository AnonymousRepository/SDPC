CXX=g++
RM=rm -f
CXXFLAGS=-std=c++11

SRCS=main.cpp DPC.cpp
OBJS=main.o DPC.o

all: SDPC clean

SDPC: $(OBJS)
	$(CXX) $(CXXFLAGS) -o SDPC $(OBJS)

DPC.o: DPC.cpp DPC.h
	$(CXX) $(CXXFLAGS) -c DPC.cpp DPC.h
main.o: main.cpp DPC.h
	$(CXX) $(CXXFLAGS) -c main.cpp DPC.h
clean:
	$(RM) $(OBJS) *.h.gch
distclean: clean
	$(RM) SDPC

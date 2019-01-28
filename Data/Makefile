CXX=g++
RM=rm -f
CXXFLAGS=-std=c++11

SRCS=processBenchmark.cpp
OBJS=processBenchmark.o

all: processBenchmark clean

processBenchmark: $(OBJS)
	$(CXX) $(CXXFLAGS) -o processBenchmark $(OBJS)

processBenchmark.o: $(SRCS)
	$(CXX) $(CXXFLAGS) -c $(SRCS)
clean:
	$(RM) $(OBJS) *.h.gch
distclean: clean
	$(RM) processBenchmark

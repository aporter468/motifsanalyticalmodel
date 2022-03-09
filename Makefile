all: main

CXXFLAGS+=-g -Wall
LDLIBS+=-lstdc++
main: Node.o Aggregator.o main.o

clean:
	rm -f *.o main


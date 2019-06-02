CC=gcc
CXX=g++
RM=rm -f


CPPFLAGS=-std=c++14 -O3
LDFLAGS=-std=c++14 -O3
LDLIBS=

PROGRAM=fress

SRCS=fress.cpp \
     main.cpp


OBJS=$(subst .cpp,.o,$(SRCS))


all: $(PROGRAM)

fress.o: fress.cpp fress.h

main.o: main.cpp fress.h


$(PROGRAM): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LDLIBS)

clean:
	rm -rf $(OBJS) $(PROGRAM)
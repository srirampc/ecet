PROG=comp2pcalign
CC=g++
CFLAGS= -c -Wall -g3 
LFLAGS=-lz
SRCS=$(wildcard *.cpp) 
OBJS=$(SRCS:.cpp=.o)

all: $(SRCS) $(PROG)

$(PROG): $(OBJS) 
	$(CC) $(OBJS) -o $@ $(LFLAGS)


$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) $<   # $< specifies the files to be compiled

clean:
	rm -rf $(PROG) $(OBJS)

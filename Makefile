# Makefile for DATCompression TopCoder Challenge

CC=g++
CFLAGS= -g
LD=g++
LDFLAGS= -g

all: compress

clean: 
	rm -f *.o
	rm compress

compress: compress.o
	$(LD) $(LDFLAGS) compress.o -o compress

compress.o: compress.cpp
	$(CC) -c $(CFLAGS) compress.cpp -o compress.o

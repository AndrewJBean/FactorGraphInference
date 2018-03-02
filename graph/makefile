CC = clang++
CFLAGS = -std=c++11 -O3

all: ImageGraph.bin

ImageGraph.bin: ImageGraph.cpp
	$(CC) $(CFLAGS) -o ImageGraph.bin -I../ ImageGraph.cpp -L./ -lfreeimage

clean:
	rm -rf *.bin

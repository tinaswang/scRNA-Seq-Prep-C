CC = g++
CFLAGS = -Wall -pedantic -ggdb -std=c++1z -g -O2
LFLAGS = -lboost_program_options
LLFLAGS = -lstdc++fs -lboost_iostreams -lpthread -lz
SRC = src
OBJ = obj
BIN = bin

barcodes: $(SRC)/get_cell_barcodes.cpp 
	$(CC) $(CFLAGS) $(SRC)/main.cpp $(LFLAGS) $(SRC)/get_cell_barcodes.hpp \
	$(SRC)/get_cell_barcodes.cpp $(LLFLAGS)

all: $(SRC)/main.cpp	
	$(SRC)/get_cell_barcodes.hpp barcodes -o preprocess
	


clean :
	rm -rf get_cell_barcodes preprocess *.dat *.o *.dSYM *.gch out

CC = g++
CFLAGS = -Wall -pedantic -ggdb -std=c++11 -g -O2
LFLAGS = -lboost_program_options
LLFLAGS = -lpthread -lz
SRC = src
OBJ = obj
BIN = bin

# barcodes: $(SRC)/get_cell_barcodes.cpp 
# 	$(CC) $(CFLAGS) $(SRC)/get_cell_barcodes.hpp \
# 	$(SRC)/get_cell_barcodes.cpp $(LLFLAGS)


# all: $(SRC)/main.cpp $(SRC)/get_cell_barcodes.hpp
# 	$(SRC)/get_cell_barcodes.hpp barcodes -o preprocess
	
all: $(SRC)/main.cpp	$(SRC)/get_cell_barcodes.hpp \
	$(SRC)/get_cell_barcodes.cpp 
	$(CC) $(CFLAGS) $(SRC)/main.cpp $(LFLAGS) $(SRC)/get_cell_barcodes.hpp \
	$(SRC)/get_cell_barcodes.cpp $(LLFLAGS) -o preprocess -lstdc++fs \
	



clean :
	rm -rf preprocess *.dat *.o *.dSYM *.gch out

/**
@file
@author Tina Wang
@date 2018
@copyright CHOOSE ONE; see License section

@brief
Functions for obtaining 

@section License

**/


#include "get_cell_barcodes.hpp"
#include <thread>
#include <set>
#include <iomanip>
vector<int> hammingCircle(string s, int n);

std::vector<int> mergeBarcodes(std::vector<int> codewords, 
						  std::vector<int> barcode_chunk);

// void writeToOutput();

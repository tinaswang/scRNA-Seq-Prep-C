/**
@file
@author Tina Wang
@date 2018
@copyright CHOOSE ONE; see License section

@brief
Functions for obtaining 

@section License

**/

// A constant value for the barcode length
#define BARCODE_LENGTH (14)
#define	UMI_LENGTH		(10)

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <string>
#include <iterator>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <future>

#include <zlib.h>
extern "C" {
#include "kseq.h"
}


std::vector<int> readBarcodes(std::string barcode_file);


int getHammingDistance(std::string s1, std::string s2);


int encode(std::string sequence); 

std::string decode(int coded);

std::unordered_map<int, int> getAllCounts(
							std::vector<std::vector<int>> barcodes);

std::vector<std::unordered_map<int, float>> getCounts(
						std::vector<std::vector<int>> barcodes);

std::vector<int> chooseBarcodesToCorrect(
									std::unordered_map<int, int> counts,
									int expected, int precomputed);

std::vector<int> getErrorCorrectionBarcodes(
								std::vector<int> codewords,
								int dmin);


void write_to_files(std::vector<std::vector<int>> barcodes, 
					
					std::vector<std::vector<int>> brc_idx_to_correct,
					int n);

int partition(std::vector<int> &barcodes, 
							std::vector<float> &values,
                           int left, int right); 

void sort(std::vector<int> &barcodes, 
				std::vector<float> &values, 
					int left, int right);

int is_far_enough(std::vector<int> codewords, 
								int i, int dmin); 

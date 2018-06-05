/**
@file file to execute our cell-barcode-getting program from
@author Tina Wang
@date 2018
@copyright Figure this out later

@brief
Functions for getting cell barcodes for the 10x pipeline. 

@section License
cool shit to fill in later

**/

#include "get_cell_barcodes.hpp"
using namespace std;

int main(int argc, char const *argv[])
{	
    vector<vector<int>> barcodes;
    vector<string> files;
    vector<unordered_map<int, float>> counts;
    vector<unordered_map<int, float>> counts_map;
	vector<vector<int>> to_correct;

	
    files.push_back("bin/testcases.fastq.gz");
    // files.push_back("../H-MEL_S1_L001_R2_001.fastq.gz");    

    // Read in all the read files
    cout << "Reading in read files" << endl;
    for (int i = 0; i < (int) files.size(); i++)
    {
    	vector<int> codes = readBarcodes(files[i]);
    	for (int j = 0; j < (int) codes.size(); j++)
    	{
    		cout << decode(codes[j]) << endl;
    	}
    	barcodes.push_back(codes);

    	cout << "read in first file" << endl;
    }

    cout << "Getting the UMI counts for each file: " << endl;
    counts = getCounts(barcodes);
    cout << "Num of unique UMIs: " << (int) counts[0].size() << endl;


    cout << "Choosing barcodes to correct" << endl;
	counts_map = chooseBarcodesToCorrect(counts, 500, 0);

    
	cout << "Getting error correction barcodes" << endl;
    to_correct = getErrorCorrectionBarcodes(counts_map, 5);
   	cout << "printing to_correct" << endl;
   	for (int i = 0; i < (int) (to_correct[0]).size(); i++)
   	{
   		cout << (to_correct[0])[i] << endl;
   	}

	int num_files = 1;
	cout << "Writing output to files" << endl;
	write_to_files(barcodes, 
					to_correct, num_files);
	return 0;
}
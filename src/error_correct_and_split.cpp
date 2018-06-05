#include "error_correct_and_split.hpp"

using namespace std;


unordered_map<int, int> neighborStrings; 

int umi_length = 10;
int barcode_length = 14;
int num_threads = 16;

unordered_map<int, int>
unordered_map<int, int> cw;
vector<int> codeword_set;
vector<int> brc_to_correct;
vector<vector<int>> barcode_split((int) num_threads, vector<int>(0));


/* @brief Generate all strings 1 char away from our input s 
 *
 * @runtime O(4n), where n is the length of s
 */
vector<int> hammingCircle(string s)
{
	vector<int> neighbors;
	alpha = 'ACTG'
	for (int i = 0; i < (int) s.length(); i++)
	{
		for (int j = 0; j < (int) alpha.length(); j++)
		{
			string updated = s;
			updated[i] = j;
			int encoded = encode(j);

			if (counts.find(encoded) != counts.end())
            {
                counts[encoded] += encode(s);
            }

            neighbors.push_back(encoded);

		}

	}

	return neighbors;
}

        

vector<int> mergeBarcodes(vector<int> barcodes, std::vector<int> codewords, 
						  std::vector<int> brc_to_correct_neigbors)
{
	int offset = barcodes[0];
	int barcs = barcodes[1];
	vector<vector<int>> retvec((int) codewords.size(), vector<int>(0));


	for(int idx = 0; idx < (int) barcodes.size(); idx++)
	{
		int barcode = barcodes[idx];
		if (codeword_set.count(barcode))
		{
			(retvec[cw[barcode]]).emplace_back(idx + offset);
		}

		else
		{
			if (brc_to_correct_neigbors.count(barcode))
			{
				neighbors = hammingCircle(decode(barcode));
				for (int i = 0; i < (int) neighbors.size(); i++)
				{
					neighbor = neighbors[i];
					if (brc_to_correct_neigbors.count(neighbor))
					{
						(retvec[cw[barcode]]).emplace_back(idx + offset);
						break;
					}
				}
			}
		}


    return retvec;
	}
}



void loadBarcodes(vector<int> barcodes, vector<int> codewords,
					vector<int> brc_idx_to_correct)

{
	set<int> s;
	unsigned int size = codewords.size();
	for(unsigned i = 0; i < size; ++i) s.insert(codewords[i]);
	codeword_set.assign( s.begin(), s.end() );

	
	vector<int> old_brc_to_correct;
	unsigned int size = brc_idx_to_correct.size();
	for(int i = 0; i < (int) size; i++)
	{
		old_brc_to_correct.emplace_back(codewords[brc_idx_to_correct[i]]);
	}
	set<int> s2;
	for(unsigned i = 0; i < size; ++i) s2.insert(old_brc_to_correct[i]);
	brc_to_correct.assign(s2.begin(), s2.end());

	
	int chunksize = (1 + (int) barcodes.size())/num_threads;

	for (int id = 0; i < (int) size; i++)
	{
		cw[codewords[id]] = id;
	}


	//  Split up the barcodes for each thread to process

	for (int i = 0; i < barcode.size(); i += chunksize)
	{
		std::vector<int> chunk;
		for (int j = i; i < i + chunksize; j++) 
   			chunk.push_back(barcodes[i]);

		(barcode_split[i]).emplace_back(chunk);
	}

	cout << "Merging barcodes" << endl;

    
  //   Generate the set of all dist-1 neighbors of brc_to_correct (for fast check in merge func)
  //   note: the number of barcodes in this set is len(brc_to_correct)*3*barcode_length

	set<int> brc_to_correct_neigbors;
	for (int i = 0; i < (int) brc_to_correct.size(); i++)
	{
		brc = brc_to_correct[i]
		vector<int> neighbors = hamming_circle(decode(brc));
		for (int j = 0; j < (int) neighbors.size(); j++) 
   			brc_to_correct_neigbors.insert(neighbors[j]);

	}

	// Create new threads
	std::thread t[num_threads];
    for (int i = 0; i < num_threads; ++i) 
    {
         t[i] = std::thread(merge_barcodes, barcode_split[i]);
    }
    cout << "Running threads" << endl;
    //Join the threads with the main thread
	for (int i = 0; i < num_threads; ++i) 
	{
		t[i].join();
   	}


   	vector<int> ret_vec;
   	for(int idd = 0; idd < (int) codewords.size(); i++)
   	{
   		// get all data returned by threads

   	}

   	vector<int> reads_per_barcode;
   	total = 0;
   	for (int i = 0; i < (int) codewords.size(); i++)
   	{
   		reads_per_barcode.push_back((int) (ret_vec[i]).size());
   		total += reads_per_barcode[i];
   	}

   	cout << "NUM_OF_READS_in_CELL_BARCODES after error-correct = " << total << endl;   
}


 /* https://stackoverflow.com/questions/667183/padding-stl-strings-in-c
  */
void padTo(std::string &str, const size_t num, const char paddingChar = '0')
{
    if(num > str.size())
        str.insert(0, num - str.size(), paddingChar);
}



void extract_all_reads(vector<string> read_files)
{

    int n_files = 30;
    int tot = (int) barcodes.size() * 4;

    double li = (double) tot / n_files;
    li = li - (li % 4) + 4;

    // You've got to read in the entire read file and split int into 30 files
	
	for (int fi = 0; fi < n_files; fi++)
	{
		int append_write = 1;
		int hi = (fi+1)*li/4;
		int lo = fi*li/4;

        string sample_name = barcode_files[0].substr(0, 
            (int) barcode_files[0].length() - end.length());


		if (fi > 0)
		{
			// change to append mode, otherwise, write new file
			append_write = 0;

			for (int cell = 0; cell < (int) codewords.size(); cell++)
			{
				string filename = sample_name + "_cell_" + padTo(std::to_string(cell), 4) + \
								"_" + decode(codewords[cell]);
				cout << "writing " << filename;
				vector<string> output_umis;
				vector<string> output_fastq;

				for (int ind = 0; ind < (int) (ret_vec[cell]).size(); ind++)
				{
					i = ret_vec[cell][ind];
					int adjusted_4i = (4 * i) % li;
					for (int l = 1; l < 5; l++)
					{
						output_fastq.push_back(all_reads at the f_i line to adjusted_4i + l)
					}
					 output_umis.push_back(ln.getline(ALL_reads_file_umi[fi],adjusted_4i+2)[BARCODE_LENGTH:BARCODE_LENGTH+UMI_LENGTH]+"\n")
				}

				// write output_umis here to filename.umi
				// write output_fastq to filename.fastq
				// record filename somewhere for the future batch file
				// compress files
			}

		}
	}
}




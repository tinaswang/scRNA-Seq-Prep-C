#include "error_correct_and_split.hpp"

using namespace std;


unordered_map<int, int> neighborStrings; 

int umi_length = 10;
int barcode_length = 14;
int num_threads = 16;

unordered_map<int, int> error_barcodes;
unordered_map<int, int> cw;
vector<int> codeword_set;
vector<int> brc_to_correct;
vector<vector<int>> barcode_split((int) num_threads, vector<int>(0));


/* @brief Generate all strings 1 char away from our input s 
 *
 * @runtime O(4n), where n is the length of s
 */
void hammingCircle(int barcode)
{
	vector<int> neighbors;
	s = decode(barcode)
	alpha = 'ACTG'
	for (int i = 0; i < (int) s.length(); i++)
	{
		for (int j = 0; j < (int) alpha.length(); j++)
		{
			string updated = s;
			updated[i] = j;
			int encoded = encode(updated);


			error_barcodes[encoded] = barcode;

		}

	}
}


        

void loadBarcodes(vector<int> &barcodes, vector<int> &codewords,
					vector<int> &brc_idx_to_correct)

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

	set<int>:: iterator it;
 	for(it = brc_to_correct.begin(); it != brc_to_correct.end(); ++it)
 	{
    	int curr = *it;
   		hammingCircle(curr);
 	}



}


vector<vector<int>> merge_barcodes(vector<int> &codewords)
{
	vector<vector<int>> retvec((int) codewords.size(), vector<int>(0));

	for (int idx = 0; idx < (int) codewords.size(); idx++)
	{

		std::vector<int>::iterator it = std::find(error_barcodes.begin(),
										 error_barcodes.end(), codewords[i]);

		// std::vector<int>::iterator it2 = std::find(brc_to_correct.begin(),
		// 								 brc_to_correct.end(), codewords[i]);

		// std::vector<int>::iterator it3 = std::find(brc_to_correct.begin(),
		// 							codeword_set(), codewords[i]);


		int barcode = codewords[idx];
		if (codeword_set.count(barcode))
		{
			(retvec[cw[barcode]]).emplace_back(idx);
		}

		else
		{
			if (it != error_barcodes.end())
			{
				neighbors = hammingCircle(decode(barcode));
				for (int i = 0; i < (int) neighbors.size(); i++)
				{
					neighbor = neighbors[i];
					if (brc_to_correct_neigbors.count(neighbor))
					{
						(retvec[cw[barcode]]).emplace_back(idx);
						break;
					}
				}
			}
		}

	}

	return retvec;
}




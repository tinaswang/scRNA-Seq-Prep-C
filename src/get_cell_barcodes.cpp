/**
@file get_cell_barcodes.cpp
@author Tina Wang
@date 2018
@copyright Figure this out later

@brief
Functions for getting cell barcodes for the 10x pipeline. 

@section License
cool shit to fill in later

**/

#include "get_cell_barcodes.hpp"
KSEQ_INIT(gzFile, gzread)

using namespace std;


/*
 * @func vector<string> readBarcodes(string barcode_file)
 *
 * @brief reads barcodes from a file
 *
 * @param string barcode_file -- the barcode file(s)? to read from
 *
 * @return vector<int> barcodes -- a vector specifying all the barcodes
 */
vector<int> readBarcodes(string barcode_file)
{
    cout << "Reading from: " << barcode_file << endl;
    vector<int> barcodes;

    gzFile fp;
    kseq_t *seq;
    int l;

    fp = gzopen(barcode_file.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {        
        string sequence = seq->seq.s;
        sequence = sequence.substr(0, BARCODE_LENGTH);
        barcodes.emplace_back(encode(sequence));
       // cout << "seq: " << sequence << endl;
    }
    // printf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp);
    return barcodes;
}




/*
 * @func getHammingDistance(string s1, string s2)
 *
 * @brief Gets hamming distance between two strings
 *  		NOTE: the two strings must be the same length!
 *
 * @param s1 -- a string
 * @param s2 -- another string
 *
 * @return dist -- an int signifying the hamming distance
 */

int getHammingDistance(string s1, string s2)
{
	int dist = 0;

 	if (s1.length() == s2.length())
    {
        for (int i=0; i< (int) s1.length(); i++)
        {
            if (s1[i] != s2[i])
            {
                dist += 1;
            }
        }
        	
    } 

    else
    {
        throw std::string("HAMMING DISTANCE ERROR: Two strings are of different length.");
    }

    return dist;
}

/*
 * @func encode(char* sequence)
 *
 * @brief Encode a DNA sequence in numerical form
 *  		
 *
 * @param sequence -- a string with the DNA sequence
 *
 * @return int coded_sequence -- a number storing the encoded sequence
 */
int encode(string sequence) 
{
	int code_sequence = 0;

	for (int i = 0; i < (int) sequence.length(); i++)
	{
		code_sequence *= 4;
		if (sequence[i] != 'N')
		{
            if (sequence[i] == 'G')
            {
			     code_sequence++;
            }
            else if (sequence[i] == 'C')
            {
                 code_sequence += 2;
            }
            else if (sequence[i] == 'T')
            {
                 code_sequence += 3;
            }
		}
		else
		{
			code_sequence += rand() % 4;
		}
	}

    return code_sequence;
}




/*
 * @func string decode(int coded)
 *
 * @brief Encode a DNA sequence in numerical form
 *  		
 * @param coded -- an int storing the encoded sequence
 *
 * @return decoded_sequence -- a string containing the decoded sequence
 *
 * 
 */
string decode(int coded)
{
	// Key for decoding letters
	string decoded_letters[4] = {"A", "G", "C", "T"}; 
	string decoded_sequence = "";
	// Iterate through and decode all letters
	for (int i = 0; i < BARCODE_LENGTH; i++)
	{
		int index = coded & 3;
		coded >>= 2;
		decoded_sequence = decoded_letters[index] + decoded_sequence;

	}
	return decoded_sequence;

}



/*
 * @func int getCounts(vector<int> barcodes)
 *
 * @brief Create an unordered_map (basically a hashtable
 *        or dictionary) of counts.
 *          
 * @param barcodes -- a vector with all the barcodes
 *
 * @return unordered_map<string, float> counts
 *
 */
vector<unordered_map<int, float>> getCounts(vector<vector<int>> barcodes)
{
    vector<unordered_map<int, float>> all_counts;


    // Iterate over the map using iterator
    for (int i = 0; i < (int) barcodes.size(); i++)
    {
        unordered_map<int, float> counts;
        vector<int> curr = barcodes[i];

        for (int j = 0; j < (int) curr.size(); j++)
        {
            if (counts.find(curr[j]) != counts.end())
            {
                counts[curr[j]] += 1;
            }   

            else 
            {
                counts[curr[j]] = 1;
            }
        }

        all_counts.push_back(counts);
    }

    return all_counts;
}

/*
 * @func int getCounts(vector<int> barcodes)
 *
 * @brief Create an unordered_map (basically a hashtable
 *        or dictionary) of counts.
 *          
 * @param barcodes -- a vector with all the barcodes
 *
 * @return unordered_map<string, float> counts
 *
 */
unordered_map<int, int> getAllCounts(vector<vector<int>> barcodes)
{
    unordered_map<int, int> counts;
    

    // Iterate over the map using iterator
    for (int i = 0; i < (int) barcodes.size(); i++)
    {
        vector<int> curr = barcodes[i];

        for (int j = 0; j < (int) curr.size(); j++)
        {
            if (counts.find(curr[j]) != counts.end())
            {
                counts[curr[j]] += 1;
            }   

            else 
            {
                counts[curr[j]] = 1;
            }
        }

    }

    return counts;
}


/*
 * @func vector<string> chooseBarcodesToCorrect(map<string, float> counts)
 *
 * @brief Remove barcodes with values outside
 *          our window parameters.
 *          
 * @param counts -- an unordered_map with barcode and count info
 *
 * @return unordered_map<string, float> counts --
 *          counts with the extraneous barcodes removed.
 *
 */

// choose barcodes for error correction
vector<int> chooseBarcodesToCorrect(
                                unordered_map<int, int> counts,
                                int expected, int precomputed)
{
    std::vector<int> labels;
    labels.reserve(counts.size());
    std::vector<float> vals;
    vals.reserve(counts.size());

    for(auto kv : counts) {
        labels.push_back(kv.first);
        vals.push_back(kv.second);  
    } 
        // std::sort(vals.begin(), vals.end(), std::greater<>());

    sort(labels, vals, 0, (int) counts.size() - 1);
    int exp_cells = (int) ((expected * 0.01) - 1);

    cout << "Num of distinct barcodes: " << (int) vals.size() << endl;

    int num_barcodes = 0;
    int num_reads_in_barcodes = 0;

    // cout << "Boundary: " << vals[exp_cells] / 10 << endl;


    // cout << "expected: " << expected << endl;
    // cout << "exp_cells: " << exp_cells << endl;

    for (int j = 0; j < (int) vals.size(); j++) 
    {
        if(vals[j] < vals[exp_cells] / 10)
        {
            num_barcodes = j;
            break;
        }

        else
        {
            num_reads_in_barcodes += vals[j];
        }
    }
    // cout << "vals[exp_cells]"<< vals[exp_cells] << endl;
    cout << "Num of barcodes: " << num_barcodes << endl;
    cout << "Num of reads in barcodes: " << num_reads_in_barcodes << endl;


    vector<int> codewords = labels;
    auto i = std::begin(codewords);
    int k = 0;
    while (i != std::end(codewords)) {
        if (k >= num_barcodes) 
            i = codewords.erase(i);
        else {
            ++i;
            k++;
            }
        }


    return codewords;
}


/* 
 * @func int is_far_enough(vector<int> barcodes, int i)
 *
 */
int is_far_enough(vector<int> codewords, int i, int dmin)
{
    vector<int> ret_vec;
    int d = dmin;

    dmin = d;
    string decoded = decode(codewords[i]);
    int j = (int) codewords.size() - 1;

    while (dmin >= d && j > 0)
    {

        if (i != j)
        {
            // compare between the hamming distance and dmin
            dmin = min(dmin, getHammingDistance(decoded, decode(codewords[j])));
        }

        j--;
    }

    if (dmin >= d)
    {
        return i;
    }

    return -1;
   
}


/*
 * @func vector<string> getErrorCorrectionBarcodes(map<string, float> counts)
 *
 * @brief return a list a barcodes to be removed
 *          
 * @param counts -- an unordered_map with barcode and count info
 *
 * @return a vector<int> of indices of barcodes to be removed
 *
 */

vector<int> getErrorCorrectionBarcodes(vector<int> codewords, int dmin)
{

    vector<int> brc_idx_to_correct;
    vector<std::future<int>> holder;
    for (int i = 0; i < (int) codewords.size(); i++)
    {
         holder.emplace_back(std::async(std::launch::async, 
                                        is_far_enough, 
                                        codewords, i, dmin));

    }

    for (auto& f : holder)
    {
        int result = (int) f.get();
        if (result != -1)
            brc_idx_to_correct.emplace_back(result);
    }

    cout << "Number of barcodes to correct: " << (int) brc_idx_to_correct.size() << endl;
    return brc_idx_to_correct;
}



/*
 * @func void write_to_files(vector<int> barcodes, 
 *                          vector<string> codewords, 
 *                          vector<int> brc_idx_to_correct)
 *
 * @brief Write vectors to output files
 *          
 * @param barcodes -- a vector with all the barcodes
 *        codewords -- a vector with all the codes
 *         brc_idx_to_correct -- a vector of integers listing
 *                               locations where we need to correct
 *
 * @return None
 *
 */

void write_to_files(vector<vector<int>> barcodes,     
					vector<vector<int>> brc_idx_to_correct, int n)
{
    // Write our barcodes to a file

    for (int i = 0; i < n; i++)
    {
        vector<int> temp_barcodes = barcodes[i];
        vector<int> temp_idx = brc_idx_to_correct[i];
        // unordered_map<int, int> temp_counts = counts[i];


        // order our counts first
        vector<int> codewords;
        codewords.reserve(temp_barcodes.size());

        for (int j = 0; j < (int) temp_barcodes.size(); j++)
        {
            codewords.push_back(temp_barcodes[j]);
        }

        if (i == 0)
        {
            ofstream output_file("barcodes.dat");
            ostream_iterator<int> output_iterator(output_file, "\n");

            copy(temp_barcodes.begin(), temp_barcodes.end(), output_iterator);

            // Write codewords to another file
            ofstream output_codes("codewords.dat");
            ostream_iterator<int> output_iterator_2(output_codes, "\n");
            copy(codewords.begin(), codewords.end(), output_iterator_2);

            // Write indexes of barcodes to correct to a file
            if ((int) brc_idx_to_correct.size() > 0)
            {
                ofstream idx_file("brc_idx_to_correct.dat");
                ostream_iterator<int> output_iterator_3(idx_file, "\n");
                copy(temp_idx.begin(), temp_idx.end(), output_iterator_3);
            }
        }

        else
        {

            ofstream b_stream("barcodes.dat", 
                                ios_base::app | ios_base::out);

            if (b_stream.fail())
            throw std::ios_base::failure(std::strerror(errno));

            //make sure write fails with exception if something is wrong
            b_stream.exceptions(b_stream.exceptions() | std::ios::failbit | std::ifstream::badbit);

            for (const auto &e : temp_barcodes) b_stream << e << "\n";

            // Write the decoded barcodes
            ofstream c_stream("codewords.dat", 
                                ios_base::app | ios_base::out);
            for (const auto &e : codewords) c_stream << e << "\n";
            if (c_stream.fail())
            throw std::ios_base::failure(std::strerror(errno));

            //make sure write fails with exception if something is wrong
            c_stream.exceptions(c_stream.exceptions() | std::ios::failbit | std::ifstream::badbit);



            // Last file
            ofstream idx("brc_idx_to_correct.dat", 
                        ios_base::app | ios_base::out);

            if (idx.fail())
            throw std::ios_base::failure(std::strerror(errno));

            //make sure write fails with exception if something is wrong
            idx.exceptions(idx.exceptions() | std::ios::failbit | std::ifstream::badbit);

            for (const auto &e : temp_idx) idx << e << "\n";

        }
    }
}


/* Helper functions */


/**
 * @brief partition() is a helper function that partitions
 *      the left and right sublists in-place for quicksort_inplace()
 *
 * @params list (the double vector), 
 *          left and right-- two ints marking the beginning and 
 *               end of the list we are currently sorting
 * @returns the highest index up to which list is sorted
 *
 *
 * algorithm: 
 *   int pivot = values[right]
 *   int i = left - 1
 *   // First, we need to rearrange the subarray in place.
 *   FOR j = left to right - 1 inclusive
 *       IF values[j] < pivot THEN
 *           i += 1
 *           swap (barcodes[i], barcodes[j])
 *       END IF
 *   END FOR
 *   IF A[hi] < A[i + 1] THEN
 *       swap (list[i + 1], list[right])
 *   END IF
 *   // Return the index for a new pivot
 *   RETURN (i + 1)
 */

int partition(vector<int> &barcodes, vector<float> &values,
                                         int left, int right) 
{
    // create the pivot
    double pivot = values[right];
    int i = left - 1;

    // rearrange the subarray in place
    for (int j = left; j < right; j++) 
    {
        if (values[j] >= pivot) 
        {
            i = i + 1;
            // if two elements are out of order, swap them
            double temp = values[i];
            values[i] = values[j];
            values[j] = temp;
            // swap the barcodes as well
            int holder = barcodes[i];
            barcodes[i] = barcodes[j];
            barcodes[j] = holder;
        }
    }
    // Check if the element at list[right] is in the wrong place
    if (values[right] >= values[i+1]) 
    {
        double temp = values[i+1];
        values[i+1] = values[right];
        values[right] = temp;
        // swap the barcodes as well
        int holder = barcodes[i+1];
        barcodes[i+1] = barcodes[right];
        barcodes[right] = holder;
    }
    return (i + 1);
}


/** 
 * @brief A modified version of quickSort used to sort the barcodes with 
 * respect to the values they make with the positive x-axis.
 * 
 * @params barcodes -- a vector of coordinates
 *          values -- a vector of doubles specifying values for each point
 *
 * @returns None; sorts in-place
 *
 * algorithm for the in-place quicksort:
 *   left: an integer indicating the leftmost index of the list we are sorting
 *   right: an integer indicating the rightmost 
 *          index of the list we are sorting
 * 
 *   IF left < right THEN
 *       int p = partition(barcodes, values, left, right)
 *       // Recursively quicksort the partitions
 *       quicksort_inplace(barcodes, values, left, p - 1 )
 *       quicksort_inplace(barcodes, values, p + 1, right)
 */

void sort(vector<int> &barcodes, vector<float> &values, 
                                        int left, int right)
{
    if (left < right) 
    {
        int p = partition(barcodes, values, left, right);
        // Sort the left partition
        sort(barcodes, values, left, p - 1);
        // Sort the right partition
        sort(barcodes, values, p + 1, right);
    }
}


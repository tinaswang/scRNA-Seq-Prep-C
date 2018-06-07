#include "get_cell_barcodes.hpp"
KSEQ_INIT(gzFile, gzread)

using namespace std;



int umi_length = 10;
// Change to 16 for v2
// Length of the cell barcode
int barcode_length = 14;
// Number of threads
int num_threads = 16;


// Keeping track of barcodes to error-correct
unordered_map<int, int> error_barcodes;

// A reverse lookup for merging barcodes
unordered_map<int, int> cw;

// The set of all unique cell barcodes
vector<int> codeword_set;

// Set of all cell barcodes to correct
vector<int> brc_to_correct;



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
        sequence = sequence.substr(0, barcode_length);
        barcodes.emplace_back(encode(sequence));

    }


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
	for (int i = 0; i < barcode_length; i++)
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
vector<unordered_map<int, float>> getCounts(vector<vector<int>> &barcodes)
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
unordered_map<int, int> getAllCounts(vector<vector<int>> &barcodes)
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
                                unordered_map<int, int> &counts,
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

    sort(labels, vals, 0, (int) counts.size() - 1);

    int automatic = vals[(counts.size()) / 2];
    int exp_cells = (int) ((expected * 0.01) - 1);
    cout << "automatic: " << automatic << endl;
    cout << "expected " << vals[exp_cells] / 10 << endl;
    cout << "Num of distinct barcodes: " << (int) vals.size() << endl;

    int num_barcodes = 0;
    int num_reads_in_barcodes = 0;


    for (int j = 0; j < (int) vals.size(); j++) 
    {
        if(vals[j] < vals[exp_cells] / 10)
        // if (vals[j] < automatic)
        {
            num_barcodes = j;
            break;
        }

        else
        {
            num_reads_in_barcodes += vals[j];
        }
    }


    cout << "Num of barcodes: " << num_barcodes << endl;
    cout << "Num of reads in barcodes: " << num_reads_in_barcodes << endl;


    vector<int> codewords = labels;
    auto i = std::begin(codewords);
    int k = 0;
    while (i != std::end(codewords)) 
    {
        if (k >= num_barcodes) 
            i = codewords.erase(i);
        else 
        {
            ++i;
            k++;
        }
    }


    return codewords;
}


/* 
 * @func int isFarEnough(vector<int> barcodes, int i)
 * 
 * @brief A helper function that does pairwise
 *          comparison between two cell barcodes
 *
 */
int isFarEnough(vector<int> &codewords, int i, int dmin = 4)
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

vector<int> getErrorCorrectionBarcodes(vector<int> &codewords, int dmin)
{

    vector<int> brc_idx_to_correct;
    vector<std::future<int>> holder;
    for (int i = 0; i < (int) codewords.size(); i++)
    {
         holder.emplace_back(std::async(std::launch::async, isFarEnough, 
                                        std::ref(codewords), i, dmin));

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




/* @brief Generate all strings 1 char away from our input s 
 *
 * @runtime O(4n), where n is the length of s
 */
void fillHammingCircleMap(int barcode)
{
    string s = decode(barcode);
    string alpha = "ACTG";
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


/* 
 * @brief Returns all strings with 1-character difference from a barcode 
 */

vector<int> hammingCircle(int barcode)
{
    vector<int> neighbors;
    string s = decode(barcode);
    string alpha = "ACTG";
    for (int i = 0; i < (int) s.length(); i++)
    {
        for (int j = 0; j < (int) alpha.length(); j++)
        {
            string updated = s;
            updated[i] = j;
            int encoded = encode(updated);
            neighbors.push_back(encoded);
        }

    }

    return neighbors;
}

/*
 * @brief Loads data into various sets
 */

        

void loadBarcodes(vector<int> &barcodes, vector<int> &codewords,
                    vector<int> &brc_idx_to_correct)

{
    set<int> s;
    int size = codewords.size();
    for(int i = 0; i < (int) size; ++i) s.insert(codewords[i]);
    codeword_set.assign( s.begin(), s.end() );

    
    vector<int> old_brc_to_correct;
    size = brc_idx_to_correct.size();
    for(int i = 0; i < (int) size; i++)
    {
        old_brc_to_correct.emplace_back(codewords[brc_idx_to_correct[i]]);
    }

    for (int id = 0; id < (int) size; id++)
    {
        cw[codewords[id]] = id;
    }


    set<int> s2;
    for(int i = 0; i < (int) size; ++i) s2.insert(old_brc_to_correct[i]);
    brc_to_correct.assign(s2.begin(), s2.end());


    for(int i = 0; i < (int) brc_to_correct.size(); i++)
    {
        fillHammingCircleMap(brc_to_correct[i]);
    }


}


/* 
 * @func mergeBarcodes(vector<int> &codewords)
 *
 * @brief Merge the actual barcodes by looking up each barcode
 *         and its neighbors
 */
vector<vector<int>> mergeBarcodes(vector<int> &codewords)
{
    int total = 0;
    vector<vector<int>> retvec((int) codewords.size(), vector<int>(0));

    for (int idx = 0; idx < (int) codewords.size(); idx++)
    {


        int barcode = codewords[idx];
        vector<int>::iterator it = find(codeword_set.begin(), 
                                            codeword_set.end(), barcode);


        if (it != codeword_set.end())
        {
            (retvec[cw[barcode]]).emplace_back(idx);
            total += 1;
        }


        else if (error_barcodes.find(codewords[idx]) != error_barcodes.end())
        {
            vector<int> neighbors = hammingCircle(barcode);
            for (int i = 0; i < (int) neighbors.size(); i++)
            {
                int neighbor = neighbors[i];
                if (error_barcodes.find(neighbor) != error_barcodes.end())
                {
                    (retvec[cw[barcode]]).emplace_back(idx);
                    std::cout << "Entered new loop" << std::endl;
                    total += 1;
                }
            }
        }
    }
        std::cout << "Total number of reads after merging: " << total << std::endl;

        return retvec;

}






/* Additional Helper functions */


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
 * 
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










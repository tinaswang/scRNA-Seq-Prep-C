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
#include <boost/program_options.hpp>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include <chrono>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
using namespace boost::program_options;
namespace po = boost::program_options;
using namespace std::chrono;

namespace fs = std::experimental::filesystem;

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} 



/*
 * has_suffix(const std::string &str, const std::string &suffix)
 * @brief Check if a string has a certain given suffix
 * @return true/false
 * https://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt
 */
bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}


int main(int argc, char **argv)
{

    try
    {
        int dmin = 5;
        int expected;
        vector<string> filepath_vec;
        // int b_length;
        // int u_length;

        // Declare the supported options.
        po::options_description desc("Allowed options: ");
        desc.add_options()
        ("help, h", "Load in the samples required")
        ("dmin, d", po::value<int>(&dmin), 
                                "set dmin for each sample like --dmin 5")
        ("expected, e", po::value<int>(&expected), "set expected value, ex. --expected 5000")
        ("filepath, f", 
            po::value<std::vector<string>>(&filepath_vec), 
            "set input directories");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        bool test = (vm.count("expected") == vm.count("filepath"));


        if (vm.count("help"))
        std::cout << desc << '\n';


        else if (test)
        {

            for (int j = 0; j < (int) filepath_vec.size(); j++)
            {
                string filepath = filepath_vec[j];
                vector<string> barcode_files;
                vector<string> read_files;


                vector<vector<int>> barcodes;
                unordered_map<int, int> counts;
                vector<int> codewords;


                DIR* dirp = opendir(filepath.c_str());
                struct dirent *dp;
                if (dirp != NULL) 
                {
                    /* print all the files and directories within directory */
                    while ((dp = readdir(dirp)) != NULL)
                    {
                        string name = dp->d_name;
                        if(has_suffix(name, "_R1_001.fastq.gz"))
                        {
                            string abs_path = filepath + name;
                            barcode_files.push_back(abs_path);
                        }


                        else if (has_suffix(dp->d_name, "_R2_001.fastq.gz"))
                        {
                            string abs_path = filepath + name;
                            read_files.push_back(abs_path);
                        }
                    }

                    closedir (dirp);
                }

                else 
                {
                    /* could not open directory */
                    perror ("");
                    return EXIT_FAILURE;
                }

                string end = "_S1_L001_R1_001.fastq.gz";
                string sample_name = barcode_files[0].substr(0, 
                    (int) barcode_files[0].length() - end.length());

                std::sort(barcode_files.begin(), barcode_files.end());
                std::sort(read_files.begin(), read_files.end());
                if ((int) barcode_files.size() != (int) read_files.size())
                {
                    cerr << "ERROR: Different number of read files and barcode files" << endl;
                }

                else if ((int) barcode_files.size() == 0)
                {
                    cerr << "ERROR: No valid files in the directory" << endl;
                }

                else
                {
                    high_resolution_clock::time_point start = high_resolution_clock::now();

                    high_resolution_clock::time_point t1 = high_resolution_clock::now();
                    cout << "Reading in barcodes" << endl;
                    for (int i = 0; i < (int) barcode_files.size(); i++)
                    {

                        vector<int> codes = readBarcodes(barcode_files[i]);
                        barcodes.push_back(codes);
                    }
                    high_resolution_clock::time_point t2 = high_resolution_clock::now();
                    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken: " << duration << endl << endl;

                    
                    cout << "Getting the raw barcode counts for each file: " << endl;
                    t1 = high_resolution_clock::now();
                    counts = getAllCounts(barcodes);
                    // int total = (int) counts.size();
                    t2 = high_resolution_clock::now();
                    duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken (microseconds): " << duration << endl << endl;

                    cout << "Choosing barcodes to correct" << endl;
                    t1 = high_resolution_clock::now();
                    codewords = chooseBarcodesToCorrect(counts, expected, 0);
                    t2 = high_resolution_clock::now();
                    duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken (microseconds): " << duration << endl << endl;

                    
                    cout << "Getting error correction barcodes w/ pairwise comparison" << endl;
                    t1 = high_resolution_clock::now();
                    vector<int> to_correct = getErrorCorrectionBarcodes(codewords, dmin);
                    t2 = high_resolution_clock::now();
                    duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken (microseconds): " << duration << endl << endl;

                    cout << "Preparing for merge..." << endl;
                    t1 = high_resolution_clock::now();
                    loadBarcodes(barcodes[j], codewords, to_correct);
                    t2 = high_resolution_clock::now();
                    duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken (microseconds): " << duration << endl << endl;

                    cout << "Merging barcodes" << endl;
                    t1 = high_resolution_clock::now();
                    vector<vector<int>> all_merged = mergeBarcodes(codewords);
                    t2 = high_resolution_clock::now();
                    duration = duration_cast<microseconds>( t2 - t1 ).count();
                    cout << "Time taken (microseconds): " << duration << endl << endl;


                    // All_merged contains all the correct barcodes -- now all
                    // we have to do is split the reads and match them up to the barcodes

                    high_resolution_clock::time_point end = high_resolution_clock::now();
                    duration = duration_cast<microseconds>(end - start).count();
                    cout << "Total Time taken (microseconds): " << duration << endl;

                }
            }
        }
    }


    catch (const error &ex)
    {
        std::cerr << ex.what() << '\n';
    }


    return 0;
}
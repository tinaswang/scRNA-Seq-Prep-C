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
#include <iostream>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;
namespace po = boost::program_options;

namespace fs = std::experimental::filesystem;

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace 
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
        vector<int> precomputed_vec;
        vector<int> dmin_vec;
        vector<int> expected_vec;
        vector<string> filepath_vec;



        // Declare the supported options.
        po::options_description desc("Allowed options: ");
        desc.add_options()
        ("help, h", "Load in the samples required")
        ("dmin, d", po::value<std::vector<int>>(&dmin_vec)->multitoken(), 
                                "set dmin for each sample like --dmin 1 2 3")
        ("precomputed, p", po::value<std::vector<int>>(&precomputed_vec), 
                                                    "set precomputed value")
        ("expected, e", po::value<std::vector<int>>(&expected_vec), "set expected value")
        ("filepath, f", 
            po::value<std::vector<string>>(&filepath_vec), 
            "set input directories")
        ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    bool test = ((vm.count("dmin") ==
                vm.count("expected")) == vm.count("filepath"));


    if (vm.count("help"))
      std::cout << desc << '\n';


    else if (vm.count("dmin") > 0 && test)
    {
        cout << "Size of dmin: " << (int) dmin_vec.size() << endl;

        cout << "Size of filepath: " << (int) filepath_vec.size() << endl;
        for (int j = 0; j < (int) filepath_vec.size(); j++)
        {
            int dmin = dmin_vec[j];
            int expected = expected_vec[j];
            string filepath = filepath_vec[j];
            // int precomputed = precomputed_vec[j];
        vector<string> barcode_files;
        vector<string> read_files;


        vector<vector<int>> barcodes;
        unordered_map<int, int> counts;
        vector<int> codewords;

        for(auto& p: fs::recursive_directory_iterator(filepath))
        {
            if(has_suffix(p.path().filename(), "_R1_001.fastq.gz"))
            {
                barcode_files.push_back(p.path());
                std::cout << p.path() << '\n';
            }

            else if (has_suffix(p.path().filename(), "_R2_001.fastq.gz"))
            {
                read_files.push_back(p.path());
                std::cout << p.path() << '\n';
            }

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

        cout << "Reading in barcodes" << endl;
        for (int i = 0; i < (int) barcode_files.size(); i++)
        {
            vector<int> codes = readBarcodes(barcode_files[i]);
            barcodes.push_back(codes);
        }



        cout << "Getting the raw barcode counts for each file: " << endl;
        counts = getAllCounts(barcodes);
        // int total = (int) counts.size();


        cout << "Choosing barcodes to correct" << endl;
        codewords = chooseBarcodesToCorrect(counts, expected, 0);
        
        cout << "Getting error correction barcodes" << endl;
        vector<int> to_correct = getErrorCorrectionBarcodes(codewords, dmin);

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
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/kseq.h>

using std::cout;
using std::endl;
using namespace std;

KSEQ_INIT(gzFile, gzread);
 
/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "get_N_runs [OPTIONS]\n");
    fprintf(stderr, "Given a FASTA file, finds the --num longest runs of Ns in each sequence and prints them in BED format.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --fasta -f The FASTA file.\n");
    fprintf(stderr, "    --num -n The number of longest N-runs to print (default: 1)\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    /** Define arguments 
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     * Fields for each argument: name, has_arg (values: no_argument, required_argument,
     *     optional_argument)
     * flag = int value to store flag for the option, or NULL if option is string
     * val = short name for string option, or NULL
     */
     
    static struct option long_options[] = {
       {"fasta", required_argument, 0, 'f'},
       {"num", required_argument, 0, 'n'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string fasta;
    int num = 1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "f:n:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'f':
                fasta = optarg;
                break;
            case 'n':
                num = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    if (fasta.length() == 0){
        fprintf(stderr, "ERROR: --fasta/-f required\n");
        exit(1);
    } 
    int progress;
    
    gzFile f_fp = gzopen(fasta.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening file %s for reading.\n", fasta.c_str());
        exit(1);
    }
   
    kseq_t* seq = kseq_init(f_fp);
    
    while ((progress = kseq_read(seq)) >= 0){
        
        // Map run length to BED-format location
        vector<pair<int, string> > runs;

        // Look for runs of N
        bool in_run = false;
        int runstart = 0;
        for (int i = 0; i < seq->seq.l; ++i){
            if (seq->seq.s[i] == 'n' || seq->seq.s[i] == 'N'){
                if (!in_run){
                    runstart = i;
                    in_run = true;
                }
            }
            else{
                if (in_run){
                    char bedstr[500];
                    sprintf(&bedstr[0], "%s\t%d\t%d", seq->name.s, runstart, i);
                    runs.push_back(make_pair(-(i - runstart), bedstr));
                    in_run = false;
                }
            }
        }
        if (in_run){
            char bedstr[500];
            sprintf(&bedstr[0], "%s\t%d\t%d", seq->name.s, runstart, seq->seq.l);
            runs.push_back(make_pair(-(seq->seq.l-runstart), bedstr));
        }
        
        sort(runs.begin(), runs.end());
        for (int i = 0; i < runs.size(); ++i){
            if (i < num){
                fprintf(stdout, "%s\n", runs[i].second.c_str());
            }
        }
    }

    // Clean up.
    kseq_destroy(seq);
    
    return 0;
}

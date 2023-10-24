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
    fprintf(stderr, "count_N_fa [OPTIONS]\n");
    fprintf(stderr, "Given a FASTA file, prints the percent N of each sequence\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --fasta -f The FASTA file.\n");
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
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string fasta;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "f:h", long_options, &option_index )) != -1){
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
        
        int ncount = 0;
        for (int i = 0; i < seq->seq.l; ++i){
            if (seq->seq.s[i] == 'n' || seq->seq.s[i] == 'N'){
                ncount++;
            }
            
        }
        fprintf(stdout, "%s\t%f\n", seq->name.s, (float)ncount/(float)seq->seq.l); 
    }

    // Clean up.
    kseq_destroy(seq);
    
    return 0;
}

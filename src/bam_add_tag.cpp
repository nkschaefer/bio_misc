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
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/sam.h>
#include <zlib.h>
#include "bam.h"

using std::cout;
using std::endl;
using namespace std;


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bam_add_tag [OPTIONS]\n");
    fprintf(stderr, "Add a tag to every entry in a BAM file. Only one tag can be given.\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --outfile -o The name of the output BAM file to create (default: stdout)\n");
    fprintf(stderr, "    --tag -t The name of the tag to create (must be 2 letters, i.e. CB)\n");
    fprintf(stderr, "    --text -T The text of the tag to add\n");
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
       {"bam", required_argument, 0, 'b'},
       {"outfile", required_argument, 0, 'o'},
       {"tag", required_argument, 0, 't'},
       {"text", required_argument, 0, 'T'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string outfile;
    string tagname;
    string tagtext;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:t:T:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'b':
                bamfile = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 't':
                tagname = optarg;
                break;
            case 'T':
                tagtext = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (bamfile.length() == 0){
        fprintf(stderr, "ERROR: bam file (--bam) required\n");
        exit(1);
    }
    if (tagname.length() != 2){
        fprintf(stderr, "ERROR: tag name must be 2 letters (i.e. CB)\n");
        exit(1);
    }
    if (tagtext.length() == 0){
        fprintf(stderr, "ERROR: tag text (-T) required\n");
        exit(1);
    }

    BGZF* outf;
    if (outfile.length() == 0){
        // Print to stdout
        outf = bgzf_open("-", "w");
    }
    else{
        outf = bgzf_open(outfile.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR opening file %s for writing.\n", outfile.c_str());
            exit(1);
        }
    }
    
    // Init BAM reader
    bam_reader reader(bamfile);
    
    // Write BAM header
    reader.write_header(outf);
    
    while(reader.next()){
        reader.add_string_tag(tagname, tagtext);
        reader.write_record(outf);
    }
    bgzf_close(outf);
    return 0;
}

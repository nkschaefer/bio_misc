/**
 * Sorts a gigantic BED file by writing out each chromosome (in sorted order) to a TMP file,
 * then merges all written TMP files.
 */
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <set>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>

using std::cout;
using std::endl;
using namespace std;

/**
 * Trim whitespace and newlines and tabs off the right end of a string.
 */
void rstrip(string& str){
    long int rpos = str.length()-1;
    while (rpos >= 0){
        if (str[rpos] == ' ' || str[rpos] == '\n' || str[rpos] == '\t'){
            rpos--;
        }
        else{
            break;
        }
    }
    if (rpos < str.length()-1){
        str = str.substr(0, rpos+1);
    }
}

/**
 * Trim whitespace and newlines and tabs off the left end of a string.
 */
void lstrip(string& str){
    long int lpos = 0;
    while (lpos < str.length()){
        if (str[lpos] == ' ' || str[lpos] == '\n' || str[lpos] == '\t'){
            lpos++;
        }
        else{
            break;
        }
    }
    if (lpos > 0){
        str = str.substr(lpos, str.length() - lpos);
    }
}

/**
 * Trim whitespace and newlines and tabs off both ends of a string.
 */
void strip(string& str){
    rstrip(str);
    lstrip(str);
}


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "sort_huge_bed [OPTIONS]\n");
    fprintf(stderr, "Sorts a gigantic BED file by writing temporary files per chromosome.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --in -i The input BED file to read\n");
    fprintf(stderr, "   --out -o The output BED file to write\n");
    fprintf(stderr, "   --tmp -t The directory to write temporary files\n");
    exit(code);
}

int main(int argc, char *argv[]) {    

    // Define arguments 
    static struct option long_options[] = {
       {"in", required_argument, 0, 'i'},
       {"out", required_argument, 0, 'o'},
       {"tmp", required_argument, 0, 't'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    
    string infile;
    string outfile;
    string tmpdir;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:t:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'i':
                infile = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 't':
                tmpdir = optarg;
                break;
            case 'h':
                help(0);
                break;
            default:
                help(0);
        }    
    }
    
    if (infile == "" || outfile == "" || tmpdir == ""){
        help(1);
    }
    
    string basename = infile;
    size_t idx = infile.rfind(".");
    if (idx != string::npos){
        basename = infile.substr(0, idx);
    }
    if (tmpdir[tmpdir.length()-1] == '/'){
        tmpdir = tmpdir.substr(0, tmpdir.length()-1);
    }
    
    // Store chromosomes in sorted order
    set<string> chroms_sorted;
    
    string prevchrom = "";
    
    // Keep current chromosome in sorted order.
    map<pair<int, int>, string> bed_curchrom;
    
    // Read input file
    fstream fin;
    fin.open(infile.c_str(), fstream::in);
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", infile.c_str());
        exit(1);
    }
    string line;
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        if (line.compare("") != 0){
            // Split by tab
            string chrom;
            int start;
            int end;
            
            istringstream tabsplitter(line);
            string fld;
            int token_index = 0;
            while(std::getline(tabsplitter, fld, '\t')){
                if (token_index == 0){
                    chrom = fld;
                }
                else if (token_index == 1){
                    start = atoi(fld.c_str());
                }
                else if (token_index == 2){
                    end = atoi(fld.c_str());
                }
                else{
                    break;
                }
                token_index++;
            }
            if (chrom != prevchrom){
                fprintf(stderr, "Processing chrom %s...\n", chrom.c_str());
                if (prevchrom != ""){
                    // Write tmp file.
                    string outname = tmpdir + "/" + basename + "_" + prevchrom + ".bed";
                    FILE* outf = fopen(outname.c_str(), "w");
                    for (map<pair<int, int>, string>::iterator bcc = bed_curchrom.begin();
                        bcc != bed_curchrom.end(); ++bcc){
                        fprintf(outf, "%s\n", bcc->second.c_str());
                    }
                    fclose(outf);
                }
                chroms_sorted.insert(chrom);
                bed_curchrom.clear();
                prevchrom = chrom;
            }
            // Insert into map.
            bed_curchrom.insert(make_pair(make_pair(start, end), line));
        }
    }
    // Handle last chromosome
    if (prevchrom != "" && bed_curchrom.size() > 0){
        // Write tmp file.
        string outname = tmpdir + "/" + basename + "_" + prevchrom + ".bed";
        FILE* outf = fopen(outname.c_str(), "w");
        for (map<pair<int, int>, string>::iterator bcc = bed_curchrom.begin();
            bcc != bed_curchrom.end(); ++bcc){
            fprintf(outf, "%s\n", bcc->second.c_str());
        }
        fclose(outf);
    }
    // Now concatenate all these sorted chunks in the correct order.
    FILE* outf = fopen(outfile.c_str(), "w");
    for (set<string>::iterator chrom = chroms_sorted.begin(); chrom != chroms_sorted.end(); ++chrom){
        string inname = tmpdir + "/" + basename + "_" + *chrom + ".bed";
        fstream inf;
        inf.open(inname.c_str(), fstream::in);
        if (!inf.good()){
            fprintf(stderr, "ERROR opening file %s\n", inname.c_str());
            exit(1);
        }
        string line;
        while (!inf.eof()){
            // Read line
            std::getline(inf, line);
            fprintf(outf, "%s\n", line.c_str());   
        }
        if (remove(inname.c_str()) != 0){
            fprintf(stderr, "ERROR deleting tmp file %s\n", inname.c_str());
        }
    }
    fclose(outf);
}

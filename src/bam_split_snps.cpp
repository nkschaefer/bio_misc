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

void parse_snps(string& fn, map<string, map<int, pair<char, char> > >& snps){
    ifstream infile(fn.c_str());
    string line;
    while(infile >> line){
        if (line.length() > 0){
            istringstream splitter(line);
            string field;
            int fld_idx = 0;
            string chrom;
            int pos;
            char allele1;
            char allele2;
            while (getline(splitter, field, '\t')){
                if (fld_idx == 0){
                    chrom = field;
                }
                else if (fld_idx == 1){
                    pos = atoi(field.c_str()) + 1;
                }
                else if (fld_idx == 3){
                    allele1 = field[0];
                }
                else if (fld_idx == 4){
                    allele2 = field[0];
                }
                ++fld_idx;
            }
            if (snps.count(chrom) == 0){
                map<int, pair<char, char> > m;
                snps.insert(make_pair(chrom, m));
            }
            snps[chrom].insert(make_pair(pos, make_pair(allele1, allele2)));
        }
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bam_split_snps [OPTIONS]\n");
    fprintf(stderr, "Takes a BAM file and a set of SNPs diagnostic of membership in one \
or another group. Outputs subsets of the BAM file containing only reads hitting one \
or the other version of each SNP. The SNP file should be in BED format (tab separated): \
chrom   pos start(0-based)  end+1(0-based)  allele1 allele2.\n");
    fprintf(stderr, "You may output BAM for both groups, or just one. One output may be \
specified as stdout (-) instead of a named output file.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --b1 -1 The output file for reads hitting allele 1 (or - if -2 is not also -)\n");
    fprintf(stderr, "    --b2 -2 The output file for reads hitting allele 2 (or - if -1 is not also -)\n");
    fprintf(stderr, "    --snps -s The BED file listing SNPs. After default BED columns, \n");
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
       {"b1", required_argument, 0, '1'},
       {"b2", required_argument, 0, '2'},
       {"snps", required_argument, 0, 's'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string out1name;
    string out2name;
    string snpfile;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:1:2:s:h", long_options, &option_index )) != -1){
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
            case '1':
                out1name = optarg;
                break;
            case '2':
                out2name = optarg;
                break;
            case 's':
                snpfile = optarg;
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
    if (out1name.length() == 0 && out2name.length() == 0){
        fprintf(stderr, "ERROR: you must provide at least one output file\n");
        exit(1);
    }
    if (out1name == "-" && out2name == "-"){
        fprintf(stderr, "ERROR: only one of -1 or -2 can be sent to stdout (-)\n");
        exit(1);
    }
    if (snpfile == ""){
        fprintf(stderr, "ERROR: SNP file is required.\n");
        exit(1);
    }
    
    map<string, map<int, pair<char, char> > > snps;
    parse_snps(snpfile, snps);
    
    // Init BAM reader
    bam_reader reader(bamfile);
    
    bool has_out1 = false;
    BGZF* out1 = NULL;
    if (out1name != ""){
        out1 = bgzf_open(out1name.c_str(), "w");
        reader.write_header(out1);
        has_out1 = true;
    }
    bool has_out2 = false;
    BGZF* out2 = NULL;
    if (out2name != ""){
        out2 = bgzf_open(out2name.c_str(), "w");
        reader.write_header(out2);
        has_out2 = true;
    }
    
    map<int, pair<char, char> >::iterator nextsnp;
    string curchrom;
    
    long int count_1 = 0;
    long int count_2 = 0;
    long int count_chimeric = 0;
    
    while(reader.next()){
        
        char* chrom = reader.ref_id();
        if (strcmp(chrom, curchrom.c_str()) != 0){
            curchrom = chrom;
            nextsnp = snps[curchrom].begin();
        }
        
        // Catch up SNPs to the current position.
        while (nextsnp != snps[curchrom].end() && nextsnp->first < reader.reference_start + 1){
            ++nextsnp;
        }
        if (nextsnp != snps[curchrom].end() && nextsnp->first >= reader.reference_start + 1 &&
            nextsnp->first <= reader.reference_end ){
            // Check all SNPs that might fall within this read.
            map<int, pair<char, char> >::iterator read_snp = nextsnp;
            // Count how many alleles the read contains that match version 1 and version 2
            // of each SNP
            int count_match1 = 0;
            int count_match2 = 0;
            while (read_snp != snps[curchrom].end() && read_snp->first >= reader.reference_start + 1 &&
                read_snp->first <= reader.reference_end){
                char allele = reader.get_base_at(read_snp->first);
                if (allele == read_snp->second.first){
                    count_match1++;
                }
                else if (allele == read_snp->second.second){
                    count_match2++;
                }
            }
            if (count_match1 > 0 && count_match2 == 0){
                count_1++;
                if (has_out1){
                    reader.write_record(out1);
                }
            }
            else if (count_match2 > 0 && count_match2 == 0){
                count_2++;
                if (has_out2){
                    reader.write_record(out2);
                }
            }
            else if (count_match1 > 0 && count_match2 > 0){
                // ???
                count_chimeric++;
            }
        }
    }
    
    fprintf(stderr, "%ld matched allele 1\n", count_1);
    fprintf(stderr, "%ld matched allele 2\n", count_2);
    fprintf(stderr, "%ld chimeric\n", count_chimeric);
    
    return 0;
}

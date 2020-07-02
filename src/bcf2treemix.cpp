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
#include <htslib/vcf.h>
#include <htslib/kseq.h>

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
 * Function to split a string by whitespace
 */
void splitstr(const string& str, vector<string>& output) { 
    istringstream buffer(str);
    copy(istream_iterator<string>(buffer), 
        istream_iterator<string>(),
        back_inserter(output));
}

char upper(char base){
    switch(base){
        case 'a':
        case 'A':
            return 'A';
            break;
        case 'c':
        case 'C':
            return 'C';
            break;
        case 'g':
        case 'G':
            return 'G';
            break;
        case 't':
        case 'T':
            return 'T';
            break;
        default:
            return base;
            break;
    }
    return base;
}

bool isbase(char base){
    if (base == 'A' || base == 'C' || base == 'G' || base == 'T'){
        return true;
    }
    return false;
}

void parse_depthfile(string& depthfile, map<string, pair<int, int> >& depths){
    ifstream infile(depthfile.c_str());
    string id;
    int mindepth;
    int maxdepth;
    while(infile >> id >> mindepth >> maxdepth){
        if (id.length() > 0){
            depths.insert(make_pair(id, make_pair(mindepth, maxdepth)));
        }
    }
}

void parse_popfile(string& popfile, map<string, string>& indv2pop){
    fstream fin;
    fin.open(popfile.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", popfile.c_str());
        exit(1);
    }
    string line;

    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        strip(line);
        if (line.length() > 0){
            istringstream tokenizer(line);
            string field;
            string indv;
            string pop;
            int token_index = 0;
            while(std::getline(tokenizer, field, '\t')){
                if (token_index == 0){
                    indv = field;
                }
                else{
                    pop = field;
                }
                token_index++;
            }
            indv2pop.insert(make_pair(indv, pop));
        }
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bcf2treemix [OPTIONS]\n");
    fprintf(stderr, "Given a BCF file (from stdin), converts it to the format \
expected by TreeMix (and prints to stdout).\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --quality -q (OPTIONAL) the minimum genotype quality required for \
a genotype to be included. If any genotype falls below this at a site, the site will be \
skipped.\n");
    fprintf(stderr, "    --depth -d (OPTIONAL) the minimum depth required for a genotype \
to be included. If any genotype falls below this at a site, the site will be excluded.\n");
    fprintf(stderr, "    --depthfile -D (OPTIONAL) A file where each line is a sample name, \
minimum coverage, and maximum coverage for that sample, space separated. This will supersede \
depth (-d) if also given.\n");
    fprintf(stderr, "    --pops -p (REQUIRED) a file mapping individual IDs to populations \
of interest. Populations missing from this file will not be included in output.\n");
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
       {"quality", required_argument, 0, 'q'},
       {"depth", required_argument, 0, 'd'},
       {"depthfile", required_argument, 0, 'D'},
       {"pops", required_argument, 0, 'p'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string depthfile;
    bool depthfile_given = false;
    int quality = -1;
    int depth = -1;
    string popfile;
    bool popfile_given = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "q:d:D:p:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'q':
                quality = atoi(optarg);
                break;
            case 'd':
                depth = atoi(optarg);
                break;
            case 'D':
                depthfile = optarg;
                depthfile_given = true;
                break;
            case 'p':
                popfile = optarg;
                popfile_given = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments
        
        map<string, pair<int, int> > depths;
        if (depthfile_given){
            parse_depthfile(depthfile, depths);
        }
        
        // Parse population file
        if (!popfile_given){
            fprintf(stderr, "ERROR: please provide a file mapping individual IDs to populations (tab separated)\n");
            exit(1);
        }
        map<string, string> indv2pop;
        parse_popfile(popfile, indv2pop);
        vector<string> indvs;
        
        // Get a set of all populations
        set<string> pops;
        for (map<string, string>::iterator ip = indv2pop.begin(); ip != indv2pop.end();
            ++ip){
            pops.insert(ip->second);
        }
        // Print in order (header)
        bool firstpop = true;
        for (set<string>::iterator p = pops.begin(); p != pops.end(); ++p){
            if (!firstpop){
                fprintf(stdout, " ");
            }
            fprintf(stdout, "%s", p->c_str());
            firstpop = false;
        }
        fprintf(stdout, "\n");
        
        string prevchrom;
        
        int num_haps = 0;
        
        int unknown_snp_id = 1;
        
        // How often should progress messages be printed?
        int progress = 50000;
        int last_printed = 0;
        
        // Convert depths to look up by index
        vector<pair<int, int> > depths_index;
        
        // Read BCF from stdin.
        bcf_hdr_t* bcf_header;
        bcf1_t* bcf_record = bcf_init();
        htsFile* bcf_reader = bcf_open("-", "r");
        if (bcf_reader == NULL){
            fprintf(stderr, "ERROR interpreting stdin as BCF format.\n");
            exit(1);
        }
        bcf_header = bcf_hdr_read(bcf_reader);
        int num_samples = bcf_hdr_nsamples(bcf_header);
        for (int i = 0; i < num_samples; ++i){
            indvs.push_back(bcf_header->samples[i]);
            if (depthfile_given){
                depths_index.push_back(depths[bcf_header->samples[i]]);
            }
        }

        fprintf(stderr, "Read %d samples\n", num_samples);
        
        // Array to store genotype data
        int32_t* gts = NULL;
        int n_gts = 0;
        
        const char* gq_key = "GQ";
        
        long int prevpos;
        
        while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
            // Make sure we're on the right chromosome.
            string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);

            if (chrom != prevchrom){
                if (last_printed > 0){
                    fprintf(stderr, "\n");
                }
                last_printed = 0;
                prevchrom = chrom;
            }
            
            // Display processing message
            if (bcf_record->pos + 1 - last_printed >= progress){
                fprintf(stderr, "Processed %s\t%d\r", chrom.c_str(), bcf_record->pos);
                last_printed = bcf_record->pos;
            }
            
            prevpos = bcf_record->pos + 1;
            
            if (bcf_record->n_allele > 2){
                // Multiallelic site; nothing to do.
                continue;
            }
            
            /*
            // Filter on variant quality?
            if (bcf_record->qual < quality){
                continue;
            }
            */
            
            // Load ref/alt alleles and other stuff
            // This puts alleles in bcf_record->d.allele[index]
            // Options for parameter 2:
            /*
            BCF_UN_STR  1       // up to ALT inclusive
            BCF_UN_FLT  2       // up to FILTER
            BCF_UN_INFO 4       // up to INFO
            BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
            BCF_UN_FMT  8                           // unpack format and each sample
            BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
            BCF_UN_ALL (BCF_UN_SHR|BCF_UN_FMT) // everything
            */
            
            bcf_unpack(bcf_record, BCF_UN_STR);
            
            bool indel = false;
            for (int i = 0; i < bcf_record->n_allele; ++i){
                if (strlen(bcf_record->d.allele[i]) > 1){
                    indel = true;
                    break;
                }
            }
            if (indel){
                //continue;
            }
            else if (!isbase(bcf_record->d.allele[0][0])){
                // Reference genome has N
                continue;
            }
            
            if (bcf_record->n_allele == 1){
                // All genotypes are homozygous reference. 
                continue;
            }
            
            // Look through all genotypes.
            
            int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
            if (num_loaded <= 0){
                fprintf(stderr, "ERROR loading genotypes at %s %d\n", 
                    chrom.c_str(), bcf_record->pos);
                exit(1);
            }
            
            set<int> gt_blacklist;
            
            // Filter out any sites where there is a genotype quality below
            // threshold (OPTIONAL)
            if (quality > 0){
                float* gqs = NULL;
                int n_gqs = 0;
                int num_gq_loaded = bcf_get_format_float(bcf_header, bcf_record, "GQ",
                    &gqs, &n_gqs);
                if (num_gq_loaded > 0){
                    for (int i = 0; i < num_samples; ++i){
                        if (!isnan(gqs[i]) && gqs[i] != bcf_float_missing &&
                            gqs[i] < quality){
                            gt_blacklist.insert(i);
                        }
                    }
                }
                free(gqs);
            }
            // Filter out any sites where there is a depth out of threshold
            // (OPTIONAL)
            if ((depth > 0 || depthfile_given)){
                int32_t* dps = NULL;
                int n_dps = 0;
                int num_dp_loaded = bcf_get_format_int32(bcf_header, bcf_record, "DP",
                    &dps, &n_dps);
                if (num_dp_loaded > 0){
                    for (int i = 0; i < num_samples; ++i){
                        if (dps[i] != bcf_int32_missing){
                            if (depthfile_given){
                                if (dps[i] < depths_index[i].first || dps[i] > depths_index[i].second){
                                    gt_blacklist.insert(i);
                                }
                            }
                            else if (dps[i] < depth){
                                gt_blacklist.insert(i);
                            }
                        }
                    }
                }
                free(dps);
            }
            
            // Store pop freqs
            map<string, int> pop2ref;
            map<string, int> pop2alt;
            for (set<string>::iterator p = pops.begin(); p != pops.end(); ++p){
                pop2ref.insert(make_pair(*p, 0));
                pop2alt.insert(make_pair(*p, 0));
            }
            
            int max_ploidy = n_gts / num_samples;
            
            for (int i = 0; i < num_samples; ++i){
                bool this_gt_missing = false;
                int32_t* gtptr = gts + i*max_ploidy;
                short count_nonref = 0;
                for (int j = 0; j < max_ploidy; ++j){
                    if (gtptr[j] == bcf_int32_vector_end){
                        // Lower ploidy
                        //geno_pass = false;
                        break;
                    }
                    else if (gt_blacklist.find(i) != gt_blacklist.end()){
                        this_gt_missing = true;
                        break;
                    }
                    else if (bcf_gt_is_missing(gtptr[j])){
                        // Missing genotype.
                        //geno_pass = false;
                        this_gt_missing = true;
                        break;
                    }
                    
                    // Don't bother looking at allele indices if all match
                    // the reference allele; we're just checking if the site
                    // should be skipped.
                    
                    if (!this_gt_missing){
                        // Retrieve allele index
                        int allele_index = bcf_gt_allele(gtptr[j]);
                        if (allele_index == 0){
                            pop2ref[indv2pop[indvs[i]]]++;
                        }
                        else{
                            pop2alt[indv2pop[indvs[i]]]++;
                        }
                    }
                }
            }
     
            // Write site to disk
            bool firstpop = true;
            for (set<string>::iterator p = pops.begin(); p != pops.end(); ++p){
                if (!firstpop){
                    fprintf(stdout, " ");
                }
                fprintf(stdout, "%d,%d", pop2ref[*p], pop2alt[*p]);
                firstpop = false;
            }
            fprintf(stdout, "\n");
            
        }
        
        bcf_hdr_destroy(bcf_header);
        bcf_destroy(bcf_record);
        bcf_close(bcf_reader);
        
        // Break carriage return line
        fprintf(stderr, "\n");
    
    return 0;
}

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

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bcf2eigenstrat [OPTIONS]\n");
    fprintf(stderr, "Given a BCF file (from stdin), converts it to the format \
expected by EIGENSTRAT and other such programs.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --output_prefix -o The base name of all output files (REQUIRED). \
Created files will have the names <prefix>.geno.gz, <prefix>.inds, <prefix>.snp, and <prefix>.rates.\n");
    fprintf(stderr, "    --ignore_phasing -p (OPTIONAL) set this flag to not check whether \
sites are phased (it will assume you know these are phased correctly).\n");
    fprintf(stderr, "    --quality -q (OPTIONAL) the minimum genotype quality required for \
a genotype to be included. If any genotype falls below this at a site, the site will be \
skipped.\n");
    fprintf(stderr, "    --depth -d (OPTIONAL) the minimum depth required for a genotype \
to be included. If any genotype falls below this at a site, the site will be excluded.\n");
    fprintf(stderr, "    --depthfile -D (OPTIONAL) A file where each line is a sample name, \
minimum coverage, and maximum coverage for that sample, space separated. This will supersede \
depth (-d) if also given.\n");
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
       {"output_prefix", required_argument, 0, 'o'},
       {"ignore_phasing", no_argument, 0, 'p'},
       {"quality", required_argument, 0, 'q'},
       {"depth", required_argument, 0, 'd'},
       {"depthfile", required_argument, 0, 'D'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string out_prefix;
    string depthfile;
    bool depthfile_given = false;
    bool ignore_phasing = false;
    int quality = -1;
    int depth = -1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:s:e:q:d:D:ph", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                out_prefix = optarg;
                break;
            case 'p':
                ignore_phasing = true;
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
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.

        // Check that output prefix was provided
        if (out_prefix.length() == 0){
            fprintf(stderr, "ERROR: you must provide a valid output prefix.\n");
            exit(1);
        }
        
        map<string, pair<int, int> > depths;
        if (depthfile_given){
            parse_depthfile(depthfile, depths);
        }
        
        // Create output files
        string out_snps_name = out_prefix + ".snp";
        string out_haps_name = out_prefix + ".ind";
        string out_geno_name = out_prefix + ".geno";
        
        FILE* out_snps_f = fopen(out_snps_name.c_str(), "w");
        if (out_snps_f == NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_snps_name.c_str());
            exit(1);
        }
        FILE* out_haps_f = fopen(out_haps_name.c_str(), "w");
        if (out_haps_f == NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_haps_name.c_str());
            exit(1);
        }
        FILE* out_geno_f = fopen(out_geno_name.c_str(), "w");
        if (!out_geno_f){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_geno_name.c_str());
            exit(1);
        }
        
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
            fprintf(out_haps_f, "%s\n", bcf_header->samples[i]);
            if (depthfile_given){
                depths_index.push_back(depths[bcf_header->samples[i]]);
            }
        }
        fclose(out_haps_f);
        fprintf(stderr, "Read %d samples\n", num_samples);
        
        // Array to store genotype data
        int32_t* gts = NULL;
        int n_gts = 0;
        
        const char* gq_key = "GQ";
        
        long int prevpos;
        
        bool ploidy_set = false;
        map<int, int> sample_ploidy;
        
        vector<short> geno_site;
        
        while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
            // Make sure we're on the right chromosome.
            string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
            
            geno_site.clear();
            
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
                continue;
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
            
            // Store which genotypes fall below quality and/or coverage thresholds
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
            if (depth > 0 || depthfile_given){
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
                        
            int max_ploidy = n_gts / num_samples;
            
            for (int i = 0; i < num_samples; ++i){
                bool this_gt_missing = false;
                int32_t* gtptr = gts + i*max_ploidy;
                short count_nonref = 0;
                for (int j = 0; j < max_ploidy; ++j){
                    if (gt_blacklist.find(i) != gt_blacklist.end()){
                        this_gt_missing = true;
                    }
                    else if (gtptr[j] == bcf_int32_vector_end){
                        // Lower ploidy
                        if (!ploidy_set && sample_ploidy.count(i) == 0){
                            // j is 1 past max ploidy index; max ploidy index is 
                            // max ploidy - 1
                            sample_ploidy.insert(make_pair(i, j));
                            if (sample_ploidy.size() == num_samples){
                                ploidy_set = true;
                            }
                        }
                        break;
                    }
                    else if (bcf_gt_is_missing(gtptr[j])){
                        // Missing genotype.
                        this_gt_missing = true;
                        break;
                    }
                    else if (!ignore_phasing && !bcf_gt_is_phased(gtptr[j])){
                        // Unphased genotype?
                        // Check whether it's homozygous (which is phased by definition).
                        bool homozygous = true;
                        for (int k = 0; k < max_ploidy; ++k){
                            if (k != j && bcf_gt_allele(gtptr[j]) != bcf_gt_allele(gtptr[k])){
                                homozygous = false;
                                break;
                            }
                        }
                        if (!homozygous){
                            gt_blacklist.insert(i);
                        }
                    }
                    
                    if (!ploidy_set && j == max_ploidy - 1  && sample_ploidy.count(i) == 0){
                        sample_ploidy.insert(make_pair(i, j+1));
                        if (sample_ploidy.size() == num_samples){
                            ploidy_set = true;
                        }
                    }
                    
                    // Don't bother looking at allele indices if all match
                    // the reference allele; we're just checking if the site
                    // should be skipped.
                    
                    // Retrieve allele index
                    int allele_index = bcf_gt_allele(gtptr[j]);
                    if (allele_index > 0){
                        count_nonref++;
                    }
                }
                
                if (this_gt_missing){
                    geno_site.push_back(9);
                }
                else{
                    if (count_nonref > 9){
                        fprintf(stderr, "ERROR: too high ploidy to fit in array (must be < 10)\n");
                        exit(1);
                    }
                    geno_site.push_back(count_nonref);
                }
            }

            // Write site to disk
            if (geno_site.size() != num_samples){
                fprintf(stderr, "ERROR: not enough sample genotypes read at site: %ld of %d\n",
                    geno_site.size(), num_samples);
                exit(1);
            }
            for (vector<short>::iterator g = geno_site.begin(); g != geno_site.end();
                ++g){
                fprintf(out_geno_f, "%d", *g);
            }
            fprintf(out_geno_f, "\n");
            
            fprintf(out_snps_f, "snp%d\t%s\t%.8f\t%d\t%c\t%c\n", unknown_snp_id,
                chrom.c_str(), (float)(bcf_record->pos + 1)*10e-8, bcf_record->pos + 1,
                bcf_record->d.allele[0][0], bcf_record->d.allele[1][0]);
            ++unknown_snp_id;
            
        }
        
        bcf_hdr_destroy(bcf_header);
        bcf_destroy(bcf_record);
        bcf_close(bcf_reader);
        free(gts);
        
        // Break carriage return line
        fprintf(stderr, "\n");
        
        // Write ploidy to file.
        string out_ploidy_name = out_prefix + ".ploidy";
        FILE* out_ploidy_f = fopen(out_ploidy_name.c_str(), "w");
        if (!out_ploidy_f){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_ploidy_name.c_str());
        }
        else{
            for (map<int, int>::iterator p = sample_ploidy.begin(); p != sample_ploidy.end();
                ++p){
                fprintf(out_ploidy_f, "%d\n", p->second);
            }
            fclose(out_ploidy_f);
        }
        
        // Clean up.
        fclose(out_snps_f);
        fclose(out_geno_f);

    
    return 0;
}

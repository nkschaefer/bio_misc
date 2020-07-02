/**
 * Computes admixture using the D/f-hat statistic from SARGE input files.
 * Uses the same parameters as admix_scan. Can be used to compare the results
 * of the two programs.
 */
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <math.h>
#include <time.h>
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
#include <zlib.h>
#include <set>

using std::cout;
using std::endl;
using namespace std;

/**
 * In the beginning, read in a chunk of the input file using the specified buffer
 * size. Find the first line in it and from it, determine the number of haplotypes
 * in the input file. Then rewind the gzFile so it can be re-read from the start.
 *
 * Having this number also tells us the number of bytes in a line in the file
 * (num_haplotypes + 1), so we can then resize our buffer to read in appropriate
 * numbers of bytes so we always end on a line break. This way, we'll never
 * have to gzseek(), which can be a very slow operation.
 */
int get_num_haplotypes(gzFile &infile, int bufsize){
    bool eof = false;
    
    while(!eof){
        char chunk[bufsize];
        int chunk_size;
        chunk_size = gzread(infile, chunk, bufsize-1);
        chunk[chunk_size] = '\0';
        eof = gzeof(infile);
        for (int i=0; i < chunk_size; i++){
            if (chunk[i] == '\n'){
                gzrewind(infile);
                return i;
            }
        }
        // Line break not found. Try again, reading in twice as much data this
        // time.
        gzrewind(infile);
        bufsize *= 2;
    }
    
    // Not found. Abort.
    fprintf(stderr, "ERROR: the number of haplotypes could not be determined from the input file.\n");
    exit(1);
}

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

int char2gt(char c){
    switch(c){
        case '0':
            return 0;
            break;
        case '1':
            return 1;
            break;
        case '2':
            return 2;
            break;
        case '3':
            return 3;
            break;
        case '4':
            return 4;
            break;
        case '5':
            return 5;
            break;
        case '6':
            return 6;
            break;
        case '7':
            return 7;
            break;
        case '8':
            return 8;
            break;
        case '9':
            return 9;
            break;
        default:
            fprintf(stderr, "ERROR: unknown genotype %c\n", c);
            exit(1);
    }
}

void process_group(string& line, set<int>& members, int& count, int& tot, vector<int>& ploidies){
    count = 0;
    tot = 0;
    for (set<int>::iterator i = members.begin(); i != members.end(); ++i){
       
        int num = char2gt(line[*i]);
        if (num != 9){
            // Not missing data.
            tot += ploidies[*i];
            count += num;
        }
    }
}

/**
 * Function to read in a chunk from the input file.
 * Returns the number of new sites added; modifies hapList, infile, and eof by reference.
 */
void read_genotypes(set<int>& ingroup,
    set<int>& unadmixed,
    set<int>& outgroup,
    set<int>& admix_in,
    set<int>& admix_out,
    vector<int>& ploidies,
    bool f_hat,
    float& admix_num_sum,
    float& admix_denom_sum,
    float& sums_numerator,
    float& sums_denominator,
    int& count_sites,
    gzFile &infile, 
    long int bufsize,
    const int num_haplotypes){
    
    long unsigned int new_sites = 0;
    char* chunk = (char*) malloc(bufsize+1);
    
    count_sites = 0;
    
    bool eof = false;
    
    int line_index = 0;
    
    while(!eof){
        int chunk_size;
        chunk_size = gzread(infile, chunk, bufsize);
        chunk[chunk_size] = '\0';
        
        eof = gzeof(infile);
        
        if (!eof){
            // Find the last newline in the chunk and rewind to just before that position.
            int steps_back = 0;
            // Should not have to do this if bufsize is properly set.
            if (chunk[chunk_size-1] != '\n'){
                for (int i=chunk_size-1; i >= 0; i--){
                    if (chunk[i] == '\n'){
                        // Truncate the string here.
                        chunk[i] = '\0';
                        gzseek(infile, -steps_back, SEEK_CUR);
                        break;
                    }
                    steps_back++;
                }
            }
            else{
                chunk[chunk_size-1] = '\0';
            }
        }
        
        // Split lines   
        istringstream linesplitter(chunk);
        string line;
        int token_index = 0;

        while(std::getline(linesplitter, line, '\n')){
            if (line.size() == num_haplotypes){
                
                new_sites++;
                
                
                int ingroup_count = 0;
                int ingroup_tot = 0;
                process_group(line, ingroup, ingroup_count, ingroup_tot, ploidies);
                if (ingroup_tot == 0){
                    ingroup_tot = 1;
                }
                
                int unadmixed_count = 0;
                int unadmixed_tot = 0;
                process_group(line, unadmixed, unadmixed_count, unadmixed_tot, ploidies);
                if (unadmixed_tot == 0){
                    unadmixed_tot = 1;
                }
                
                int admixer_in_count = 0;
                int admixer_in_tot = 0;
                process_group(line, admix_in, admixer_in_count, admixer_in_tot, ploidies);
               
                int admixer_out_count = 0;
                int admixer_out_tot = 0;
                process_group(line, admix_out, admixer_out_count, admixer_out_tot, ploidies);
                
                int admixer_all_tot = admixer_in_tot + admixer_out_tot;
                if (admixer_all_tot == 0){
                    admixer_all_tot = 1;
                }
                if (admixer_in_tot == 0){
                    admixer_in_tot = 1;
                }
                if (admixer_out_tot == 0){
                    admixer_out_tot = 1;
                }
                
                int outgroup_count = 0;
                int outgroup_tot = 0;
                process_group(line, outgroup, outgroup_count, outgroup_tot, ploidies);
                if (outgroup_tot == 0){
                    outgroup_tot = 1;
                }
                
                float p1 = (float)unadmixed_count / (float)unadmixed_tot;
                float p2 = (float)ingroup_count / (float)ingroup_tot;
                float p3 = (float)(admixer_in_count + admixer_out_count) / (float)(admixer_all_tot);
                float p4 = (float)outgroup_count / (float)outgroup_tot;
                
                sums_numerator += ((1-p1)*p2*p3*(1-p4) - p1*(1-p2)*p3*(1-p4));
                sums_denominator += ((1-p1)*p2*p3*(1-p4) + p1*(1-p2)*p3*(1-p4));
                
                if (f_hat){
                    float p3_1 = (float)admixer_in_count / (float)admixer_in_tot;
                    float p3_2 = (float)admixer_out_count / (float)admixer_out_tot;
                    admix_num_sum += ((1-p1)*p3_1*p3_2*(1-p4) - p1*(1-p3_1)*p3_2*(1-p4));
                    admix_denom_sum += ((1-p1)*p3_1*p3_2*(1-p4) + p1*(1-p3_1)*p3_2*(1-p4));
                }                
                count_sites++;
                line_index++;
            }
        }
        // Avoid letting garbage into this buffer next time
        memset(chunk, 0, bufsize);
    }
    free(chunk);
}

/**
 * Parse the list of individuals.
 **/
void parse_indvs(vector<string>& indvs, string filename){
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        if (line.compare("") != 0){
            // Store in vector
            indvs.push_back(line);
        }
    }
}

/**
 * Parse the list of individuals, but return a mapping of name -> haplotype index,
 * rather than the other way around
 */
void parse_indvs_map(map<string, int>& indvs, string filename){
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    
    int hap_index = 0;
    
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        if (line.compare("") != 0){
            indvs.insert(make_pair(line, hap_index));
            hap_index++;
        }
    }
}

/**
 * Parses a file mapping individual to population, tab-separated.
 */
void parse_pops(map<string, string>& indv2pop, map<string, vector<string> >& pop2indv, string filename){
    fstream fin;
    fin.open(filename.c_str(), fstream::in);
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    while (!fin.eof()){
        std::getline(fin, line);
        if (line.compare("") != 0){
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
            
            string indvname;
            string popname;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    indvname = field;
                }
                else if (field_index == 1){
                    popname = field;
                }
                field_index++;
            }
            indv2pop.insert(make_pair(indvname, popname));
            if (pop2indv.count(popname) == 0){
                vector<string> indvs;
                pop2indv.insert(make_pair(popname, indvs));
            }
            pop2indv[popname].push_back(indvname);          
        }
    }
}


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "eig_dstat [OPTIONS]\n");
   fprintf(stderr, "Computes admixture (using same parameters as admix_scan) \
but uses the D/f-hat statistics on an eigenstrat geno file.\n");
   fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --geno -g Input genotype file (can be gzipped)\n");
   fprintf(stderr, "    --ingroup -i One or more populations that are admixed (can \
specify more than once)\n");
    fprintf(stderr, "   --nonadmixed -n One or more populations related to the ingroup, but \
believed to be unadmixed (can specify more than once)\n");
   fprintf(stderr, "    --outgroup -o One or more populations that are not admixed (can \
specify more than once)\n");
   fprintf(stderr, "    --admixer -a The admixing population of interest\n");
fprintf(stderr, "    --indvs -v The file containing names of individuals \n");
    fprintf(stderr, "    --pops -p A file mapping individual names to population names, \
tab-separated\n");
    fprintf(stderr, "    --proportion -f Specify to calculate admixture proportion per \
population instead of D statistic. Will randomly hold out 50%% of outgroup haplotypes\n");
    fprintf(stderr, "    --ploidy -P The file listing ploidy of each sample");
    fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
exit(code);
}

// Define a structure to store summary data about admixture within a population.
struct pop_stats {
    string popname;
    float perc_admixed;
    float perc_admixed_sites;
    float tract_len;
    long int tract_len_min;
    long int tract_len_max;
    float allele_freq;
};

bool pop_stats_comp(pop_stats s1, pop_stats s2){
    if (s1.perc_admixed == s2.perc_admixed){
        return (s1.popname.compare(s2.popname) > 0);
    }
    else{
        return s1.perc_admixed > s2.perc_admixed;
    }
}


void parse_ploidies(string& filename, vector<int>& ploidy){
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;

    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        strip(line);
        if (line.length() > 0){
            int p = atoi(line.c_str());
            ploidy.push_back(p);
        }
    }
}

int main(int argc, char *argv[]) {    
    // Define arguments 
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
    // Fields for each argument: name, has_arg (values: no_argument, required_argument,
    //     optional_argument)
    // flag = int value to store flag for the option, or NULL if option is string
    // val = short name for string option, or NULL
    static struct option long_options[] = {
       {"ingroup", required_argument, 0, 'i'},
       {"nonadmixed", required_argument, 0, 'n'},
       {"outgroup", required_argument, 0, 'o'},
       {"admixer", required_argument, 0, 'a'},
       {"indvs", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"bufsize", optional_argument, 0, 'b'},
       {"geno", required_argument, 0, 'g'},
       {"proportion", no_argument, 0, 'f'},
       {"ploidy", required_argument, 0, 'P'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    
    vector<string> ingroup_pops;
    vector<string> nonadmixed_pops;
    vector<string> outgroup_pops;
    string admixer;
    string indvfilename;
    string popfilename;
    string genofilename;
    string ploidyfile;
    bool ploidy_given = false;
    bool f_hat = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:a:v:p:b:g:n:P:fh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'i':
                ingroup_pops.push_back(string(optarg));
                break;
            case 'n':
                nonadmixed_pops.push_back(string(optarg));
                break;
            case 'o':
                outgroup_pops.push_back(string(optarg));
                break;
            case 'a':
                admixer = string(optarg);
                break;
            case 'v':
                indvfilename = string(optarg);
                break;
            case 'p':
                popfilename = string(optarg);
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'g':
                genofilename = optarg;
                break;
            case 'P':
                ploidyfile = optarg;
                ploidy_given = true;
                break;
            case 'f':
                f_hat = true;
                break;
            case '?':
                //help(0);
                break;
            case 'h':
                help(0);
                break;
            default:
                help(0);
        }    
    }
    
    if (ingroup_pops.size() == 0 || outgroup_pops.size() == 0 || admixer.length() == 0){
        fprintf(stderr, "ERROR: you must provide at least one ingroup, outgroup, and admixing population.\n");
        exit(1);
    }
    if (genofilename.length() == 0){
        fprintf(stderr, "ERROR: input genotype file name not provided.\n");
        exit(1);
    }
    if (!ploidy_given){
        fprintf(stderr, "ERROR: missing ploidy input file\n");
        exit(1);
    }
    
    vector<int> ploidies;
    parse_ploidies(ploidyfile, ploidies);
    
    // Get file reading stuff ready
    FILE *instream = stdin;
    if (instream == NULL){
        fprintf(stderr, "Error opening input stream\n");
        exit(1);
    }
    gzFile fp = gzdopen(fileno(instream), "rb");
    if (!fp){
        fprintf(stderr, "ERROR: unable to read from stdin.\n");
        exit(1);
    }
    
    vector<string> indvs;
    parse_indvs(indvs, indvfilename);
    
    map<string, int> indv2hap;
    parse_indvs_map(indv2hap, indvfilename);
    map<string, string> indv2pop;
    map<string, vector<string> > pop2indv;
    parse_pops(indv2pop, pop2indv, popfilename);
    
    vector<string> leafname_pops;
    for (vector<string>::iterator indvit = indvs.begin(); indvit != indvs.end();
        ++indvit){
        //string indv_sub = (*indvit).substr(0, (*indvit).find('-'));
        leafname_pops.push_back(indv2pop[*indvit]);
    }
    
    gzFile geno_in = gzopen(genofilename.c_str(), "r");
    if (!geno_in){
        fprintf(stderr, "ERROR opening %s for reading.\n", genofilename.c_str());
        exit(1);
    }
    // Ensure that the buffer is big enough to read in at least one whole row of 
    // haplotype data from the input file at a time
    int num_haplotypes = get_num_haplotypes(geno_in, bufsize);
    if (num_haplotypes != indvs.size()){
        fprintf(stderr, "Differing number of samples in input geno file and sample list.\n");
        exit(1);
    }
    
    fprintf(stderr, "%d samples in file\n", num_haplotypes);
    fprintf(stderr, "%ld sample names read\n", indvs.size());
    
    fprintf(stderr, "D(");
    for (int i = 0; i < nonadmixed_pops.size(); ++i){
        fprintf(stderr, "%s", nonadmixed_pops[i].c_str());
    }
    fprintf(stderr, ", ");
    for (int i = 0; i < ingroup_pops.size(); ++i){
        fprintf(stderr, "%s", ingroup_pops[i].c_str());
    }
    fprintf(stderr, ", ");
    fprintf(stderr, "%s", admixer.c_str());
    fprintf(stderr, ", ");
    for (int i = 0; i < outgroup_pops.size(); ++i){
        fprintf(stderr, "%s", outgroup_pops[i].c_str());
    }
    fprintf(stderr, ")\n");
    
    // We'll now make sure that bufsize is a multiple of num_haplotypes+1 so 
    // we never have to gzseek().
    if (bufsize < num_haplotypes + 1){
        bufsize = num_haplotypes + 1;
    }
    else{
        bufsize = (long int) round( (float) (num_haplotypes+1) * round((float) bufsize / (float) (num_haplotypes+1)));
    }
    
    // Initialize random number seed    
    srand(time(NULL));
    
    float sums_numerator = 0.0;
    float sums_denominator = 0.0;
    // Store numerator sum for outgroup subsample, so we can calculate f hat
    float admix_num_sum = 0.0;
    float admix_denom_sum = 0.0;
    
    long int count_ig_haps = 0;
    
    set<int> ingroup_set;
    for (vector<string>::iterator ingroup_pop = ingroup_pops.begin();
        ingroup_pop != ingroup_pops.end(); ++ingroup_pop){
        for (vector<string>::iterator ingroup_name = pop2indv[*ingroup_pop].begin();
            ingroup_name != pop2indv[*ingroup_pop].end(); ++ingroup_name){
            if (indv2hap.count(*ingroup_name) > 0){
                ingroup_set.insert(indv2hap[*ingroup_name]);
            }
        }
    }
    
    set<int> unadmixed_set;
    for (vector<string>::iterator unadm_pop = nonadmixed_pops.begin();
        unadm_pop != nonadmixed_pops.end(); ++unadm_pop){
        for (vector<string>::iterator unadm_name = pop2indv[*unadm_pop].begin();
            unadm_name != pop2indv[*unadm_pop].end(); ++unadm_name){
            if (indv2hap.count(*unadm_name) > 0){
                unadmixed_set.insert(indv2hap[*unadm_name]);
            }
        }
    }
    
    set<int> outgroup_set;
    for (vector<string>::iterator outgroup_pop = outgroup_pops.begin();
        outgroup_pop != outgroup_pops.end(); ++outgroup_pop){
        for (vector<string>::iterator outgroup_name = pop2indv[*outgroup_pop].begin();
            outgroup_name != pop2indv[*outgroup_pop].end(); ++outgroup_name){
            if (indv2hap.count(*outgroup_name) > 0){
                outgroup_set.insert(indv2hap[*outgroup_name]);
                //fprintf(stderr, "outgroup %d %s\n", indv2hap[*outgroup_name], outgroup_name->c_str());
            }
        }
    }

    set<int> admix_haps_all;
    for (vector<string>::iterator admix_hap = pop2indv[admixer].begin();
        admix_hap != pop2indv[admixer].end(); ++admix_hap){
        if (indv2hap.count(*admix_hap) > 0){
            admix_haps_all.insert(indv2hap[*admix_hap]);
        }
    }
    
    // Clades for calculating f hat
    set<int> admix_in_set;
    set<int> admix_out_set;
   
    if (f_hat){
        // Sort randomly
        vector<int> admix_unsorted;
        for (set<int>::iterator a = admix_haps_all.begin(); a != admix_haps_all.end();
            ++a){
            admix_unsorted.push_back(*a);
        }
        random_shuffle(admix_unsorted.begin(), admix_unsorted.end());
        for (int i = 0; i < (int)floor((float)admix_haps_all.size() / 2.0); ++i){
            admix_in_set.insert(admix_unsorted[i]);
        }
        for (int i = admix_in_set.size(); i < admix_unsorted.size(); ++i){
            admix_out_set.insert(admix_unsorted[i]);
        }
        
        if (admix_in_set.size() < 1 || admix_out_set.size() < 1){
            fprintf(stderr, "ERROR: unable to split admixers into two subpopulations for \
calculating f-hat. Please re-run with more haplotypes.\n");
            exit(1);
        }
       
        fprintf(stderr, "Split admixers into %ld ingroup and %ld outgroup samples for f-hat calculation\n",
            admix_in_set.size(), admix_out_set.size());
    }
    else{
        admix_in_set = admix_haps_all;
    }
    
    if (ingroup_set.size() == 0){
        fprintf(stderr, "ERROR: no ingroup individuals found\n");
        exit(1);
    }
    if (unadmixed_set.size() == 0){
        fprintf(stderr, "ERROR: no nonadmixed individuals found\n");
        exit(1);
    }
    if (admix_in_set.size() == 0 && admix_out_set.size() == 0){
        fprintf(stderr, "ERROR: no admixing individuals found\n");
        exit(1);
    }
    if (outgroup_set.size() == 0){
        fprintf(stderr, "ERROR: no outgroup individuals found\n");
        exit(1);
    }
    
    int num_sites = 0;

    read_genotypes(ingroup_set, unadmixed_set, outgroup_set, admix_in_set, admix_out_set, ploidies, 
        f_hat, admix_num_sum, admix_denom_sum, sums_numerator, sums_denominator, num_sites, geno_in, bufsize, num_haplotypes);
    gzclose(geno_in);
    
    // Print results.
    fprintf(stdout, "== Population-level D-statistics ==\n");
   
    float d = sums_numerator / sums_denominator;
    if (f_hat){
        float f = (sums_numerator / sums_denominator) / (admix_num_sum / admix_denom_sum);
        fprintf(stdout, "%f\t%f\n", d, f);
    }
    else{
        fprintf(stdout, "%f\n", d);
    }
    
}

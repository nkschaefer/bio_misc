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

void parse_snps(string& fn, map<string, map<int, char> >& snps1,
    map<string, map<int, char> >& snps2, bool fixed_only){
    ifstream infile(fn.c_str());
    
    string chrom;
    string pos;
    string pos2;
    string alleles1;
    string alleles2;
    
    long int snpcount = 0;
    
    long int snps1count = 0;
    long int snps2count = 0;

    while(infile >> chrom >> pos >> pos2 >> alleles1 >> alleles2){
        if (fixed_only){
            if (alleles1.length() == 1 && alleles2.length() == 1 && alleles1[0] != alleles2[0]){
                if (snps1.count(chrom) == 0){
                    map<int, char> m;
                    snps1.insert(make_pair(chrom, m));
                }
                snps1[chrom].insert(make_pair(atoi(pos.c_str())+1, alleles1[0]));
                snps1count++;
                if (snps2.count(chrom) == 0){
                    map<int, char> m;
                    snps2.insert(make_pair(chrom, m));
                }
                snps2[chrom].insert(make_pair(atoi(pos.c_str())+1, alleles2[0]));
                snps2count++;              
            }
        }
        else{
            for (int i = 0; i < alleles1.length(); ++i){
                if (alleles1[i] != ','){
                    bool uniq = true;
                    for (int j = 0; j < alleles2.length(); ++j){
                        if (alleles2[j] != ',' && alleles2[j] == alleles1[i]){
                            uniq = false;
                            break;
                        }
                    }
                    if (uniq){
                        // Diagnostic of species 1
                        if (snps1.count(chrom) == 0){
                            map<int, char> m;
                            snps1.insert(make_pair(chrom, m));
                        }
                        snps1[chrom].insert(make_pair(atoi(pos.c_str())+1, alleles1[i]));
                        snps1count++;
                    }
                }
            }
            for (int i = 0; i < alleles2.length(); ++i){
                if (alleles2[i] != ','){
                    bool uniq = true;
                    for (int j = 0; j < alleles1.length(); ++j){
                        if (alleles1[j] != ',' && alleles1[j] == alleles2[i]){
                            uniq = false;
                            break;
                        }
                    }
                    if (uniq){
                        // Diagnostic of species 2
                        if (snps2.count(chrom) == 0){
                            map<int, char> m;
                            snps2.insert(make_pair(chrom, m));
                        }
                        snps2[chrom].insert(make_pair(atoi(pos.c_str())+1, alleles2[i]));
                        snps2count++;
                    }
                }
            }
        }
        ++snpcount;
        if (snpcount % 10000 == 0){
            fprintf(stderr, "%ld SNPs read\r", snpcount);
        }
    }
    fprintf(stderr, "\nSpecies 1: %ld informative SNPS | Species 2: %ld informative SNPs\n", snps1count, snps2count);
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
    fprintf(stderr, "    --paired -p Specify paired reads. Here, if either read in a pair is species-specific, \
then both reads will be taken as mapping to that species. This requires 2 passes of the file.\n");
    fprintf(stderr, "    --alleleHMM -a Instead of splitting into BAM files, output alelle counts \
in the format expected by the AlleleHMM program\n");
    fprintf(stderr, "    --cell_barcode -c If using 10x Cell Ranger output, specify the output data file \
listing valid barcodes (one per line). This program will read all possible barcodes from that file and will add the cell barcode \
(CB tag) to each depth count for each parental allele of each SNP. Valid only if using alleleHMM output.\n");
    fprintf(stderr, "    --window -w Rather than splitting BAM or making AlleleHMM output, output number of reads \
specific to each species within genomic windows of this size. Overrides -a, -1, and -2.\n");
    fprintf(stderr, "    --offset -o If using windows, this is how far offset each new window should be (omitted = \
non-overlapping windows)\n");
    fprintf(stderr, "    --bed -B Instead of using fixed-width windows, provide a BED file. Counts per cell of \
reads overlapping each feature in the BED file, hitting species 1 alleles and hitting species 2 alleles will be \
output\n");
    fprintf(stderr, "    --fixed -f only consider SNPs fixed in the two populations (i.e. neither allele \
is shared)\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

void print_alleleHMM(string& curchrom, 
    unordered_map<string, map<int, int> >& counts1, 
    unordered_map<string, map<int, int> >& counts2, 
    set<int>& allkeys, 
    int curpos,
    bool cell_barcode){
    for (set<int>::iterator k = allkeys.begin(); k != allkeys.end(); ){
        if (curpos == -1 || *k < curpos){
            
            for (unordered_map<string, map<int, int> >::iterator it1 = counts1.begin();
                it1 != counts1.end(); ++it1){
                int c1 = 0;
                if (it1->second.count(*k) > 0){
                    c1 = it1->second[*k];
                    it1->second.erase(*k);
                } 
                int c2 = 0;
                if (counts2[it1->first].count(*k) > 0){
                    c2 = counts2[it1->first][*k];
                    counts2[it1->first].erase(*k);
                }
                // Avoid printing lines with no counts, unless this is bulk data
                if (!cell_barcode || c1 > 0 || c2 > 0){
                    fprintf(stdout, "%s\t%d\t%d\t%d", curchrom.c_str(), *k, c1, c2);
                    if (cell_barcode){
                        fprintf(stdout, "\t%s", it1->first.c_str());
                    }
                    fprintf(stdout, "\n");
                }
            }
            allkeys.erase(k++);
        }
        else{
            break;
        }
    }
}

void parse_bed_file(string& filename, map<string, set<pair<int, int> > >& bed){
    ifstream infile(filename);
    string line;
    while (getline(infile, line)){
        istringstream splitter(line);
        string field;
        string chrom;
        int start;
        int end;
        int idx = 0;
        while (getline(splitter, field, '\t')){
            if (idx == 0){
                chrom = field;
            }
            else if (idx == 1){
                start = atoi(field.c_str());
            }
            else if (idx == 2){
                end = atoi(field.c_str());
            }
            else{
                break;
            }
            ++idx;
        }
        if (bed.count(chrom) == 0){
            set<pair<int, int> > s;
            bed.insert(make_pair(chrom, s));
        }
        bed[chrom].insert(make_pair(start+1, end));
    }   
}

void parse_barcode_file(string& filename, set<string>& all_barcodes){
    ifstream infile(filename);
    bool firstline = true;
    string line;
    while (getline(infile, line)){
        if (firstline){
            firstline = false;
        }
        else{
            /*
            istringstream splitter(line);
            string field;
            string cellid;
            while(getline(splitter, field, ',')){
                cellid = field;
                break;
            }
            if (cellid != "NO_BARCODE"){
                all_barcodes.insert(cellid);
            }
            */
            all_barcodes.insert(line);
        }   
    }
}

/**
 *  Print out data from windows currently being tracked and remove what's possible from
 *  data structures. A position of -1 means the entire chromosome is finished.
 */
void process_windows(unordered_map<string, map<int, int> >& counts1, 
    unordered_map<string, map<int, int> >& counts2,
    set<string>& all_barcodes,
    string& chrom, 
    long int pos, 
    int window, 
    int offset, 
    set<int>& winstarts, 
    bool has_barcodes){
    for (set<int>::iterator start = winstarts.begin(); start != winstarts.end(); ){
        if (pos == -1 || pos > *start + window){
            for (set<string>::iterator bc = all_barcodes.begin(); bc != all_barcodes.end();
                ++bc){
                int count1 = 0;
                if (counts1.count(*bc) > 0){
                    if (counts1[*bc].count(*start) > 0){
                        count1 = counts1[*bc][*start];
                        counts1[*bc].erase(*start);
                    }
                }
                int count2 = 0;
                if (counts2.count(*bc) > 0){
                    if (counts2[*bc].count(*start) > 0){
                        count2 = counts2[*bc][*start];
                        counts2[*bc].erase(*start);
                    }
                }
                if (count1 > 0 || count2 > 0){
                    if (has_barcodes){
                        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%s\n", chrom.c_str(), 
                            *start, *start + window, count1, count2, bc->c_str());
                    }
                    else{
                        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\n", chrom.c_str(),
                            *start, *start + window, count1, count2);
                    }
                }
            }
            winstarts.erase(start++);
        }
        else{
            break;
        }
    }   
}

void process_bed(unordered_map<string, map<pair<int, int>, int> >& bedcounts1,
    unordered_map<string, map<pair<int, int>, int> >& bedcounts2,
    set<string>& all_barcodes,
    string& chrom,
    int pos,
    map<string, set<pair<int, int> > >& intervals,
    bool cell_barcode){
    
    if (intervals.count(chrom) == 0){
        return;
    }

    for (set<pair<int, int> >::iterator curbed = intervals[chrom].begin(); curbed != intervals[chrom].end(); ){
        if (pos == -1 || curbed->second < pos){
            // Done with this.
            for (set<string>::iterator bc = all_barcodes.begin(); bc != all_barcodes.end(); ++bc){
                int count1 = 0;
                int count2 = 0;
                if (bedcounts1.count(*bc) > 0){
                    if (bedcounts1[*bc].count(*curbed) > 0){
                        count1 = bedcounts1[*bc][*curbed];
                    }
                    bedcounts1[*bc].erase(*curbed);
                }
                if (bedcounts2.count(*bc) > 0){
                    if (bedcounts2[*bc].count(*curbed) > 0){
                        count2 = bedcounts2[*bc][*curbed];
                    }
                    bedcounts2[*bc].erase(*curbed);
                }
                if (count1 > 0 || count2 > 0){
                    if (!cell_barcode){
                        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\n", chrom.c_str(),
                            curbed->first-1, curbed->second, count1, count2);
                    }
                    else{
                        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%s\n", chrom.c_str(),
                            curbed->first-1, curbed->second, count1, count2, bc->c_str());
                    } 
                }
            }       
            intervals[chrom].erase(curbed++);   
        }
        else if (pos != -1 && curbed->first > pos){
            // In the future.
            break;
        }
        else{
            ++curbed;
        }
    }

    if (pos == -1){
        // Clear everything.
        intervals.erase(chrom);   
    }   
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
       {"paired", no_argument, 0, 'p'},
       {"alleleHMM", no_argument, 0, 'a'},
       {"cell_barcode", required_argument, 0, 'c'},
       {"window", required_argument, 0, 'w'},
       {"offset", required_argument, 0, 'o'},
       {"bed", required_argument, 0, 'B'},
       {"fixed", no_argument, 0, 'f'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string out1name;
    string out2name;
    string snpfile;
    bool alleleHMM = false;
    bool paired = false;
    bool cell_barcode = false;
    string cell_barcode_file = "";
    int window = -1;
    int offset = -1;
    bool fixed = false;
    bool bedfile_given = false;
    string bedfile;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:1:2:s:c:w:o:B:fpah", long_options, &option_index )) != -1){
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
            case 'p':
                paired = true;
                break;
            case 'B':
                bedfile_given = true;
                bedfile = optarg;
                break;
            case 'f':
                fixed = true;
                break;
            case 'a':
                alleleHMM = true;
                break;
            case 'c':
                cell_barcode = true;
                cell_barcode_file = optarg;
                break;
            case 'w':
                window = atoi(optarg);
                break;
            case 'o':
                offset = atoi(optarg);
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
    if (!alleleHMM && window == -1 && !bedfile_given){
        if (out1name.length() == 0 && out2name.length() == 0){
            fprintf(stderr, "ERROR: you must provide at least one output file\n");
            exit(1);
        }
        if (out1name == "-" && out2name == "-"){
            fprintf(stderr, "ERROR: only one of -1 or -2 can be sent to stdout (-)\n");
            exit(1);
        }
    }
    if (cell_barcode && !alleleHMM && window == -1 && !bedfile_given){
        fprintf(stderr, "ERROR: cell_barcode argument only applicable with alleleHMM, window, or BED output\n");
        exit(1);
    }
    if (snpfile == ""){
        fprintf(stderr, "ERROR: SNP file is required.\n");
        exit(1);
    }
    short opts_chosen = 0;
    if (alleleHMM){
        opts_chosen++;
    }
    if (window != -1){
        opts_chosen++;
    }
    if (bedfile_given){
        opts_chosen++;
    }

    if (opts_chosen > 1){
        fprintf(stderr, "ERROR: only one of alleleHMM, window, or BED output is allowed\n");
        exit(1);
    } 
    if (window != -1 && offset != -1 && offset > window){
        fprintf(stderr, "ERROR: offset can not be greater than window size.\n");
        exit(1);
    }
    if (window != -1 && offset <= 0){
        offset = window;
    }
    map<string, set<pair<int, int> > > bed;
    if (bedfile_given){
        parse_bed_file(bedfile, bed);
    }
    set<string> all_barcodes;
    if (cell_barcode){
        parse_barcode_file(cell_barcode_file, all_barcodes);
    }
    else{
        all_barcodes.insert("");
    }
    map<string, map<int, char> > snps1;
    map<string, map<int, char> > snps2;
    parse_snps(snpfile, snps1, snps2, fixed);
    
    // To track paired reads
    unordered_map<string, int> read_count_snp1;
    unordered_map<string, int> read_count_snp2;

    bool has_out1 = false;
    BGZF* out1 = NULL;
    if (out1name != ""){
        out1 = bgzf_open(out1name.c_str(), "w");
        has_out1 = true;
    }
    bool has_out2 = false;
    BGZF* out2 = NULL;
    if (out2name != ""){
        out2 = bgzf_open(out2name.c_str(), "w");
        has_out2 = true;
    }
    
    long int tot = 0;
    long int count_1 = 0;
    long int count_2 = 0;
    long int count_chimeric = 0;
    long int count_none = 0;
    
    if (!alleleHMM){
        paired = true;
    } 
    // Store counts (if alleleHMM format - mapped to barcode, or empty string if no barcodes)
    // If counting in windows, each first int key will be a start coordinate of a window, rather
    // than a SNP position
    unordered_map<string, map<int, int> > counts1;
    unordered_map<string, map<int, int> > counts2;
    unordered_map<string, map<pair<int, int>, int> > bedcounts1;
    unordered_map<string, map<pair<int, int>, int> > bedcounts2;
    for (set<string>::iterator bc = all_barcodes.begin(); bc != all_barcodes.end(); ++bc){
        map<int, int> m;
        counts1.insert(make_pair(*bc, m));
        counts2.insert(make_pair(*bc, m));
        if (bedfile_given){
            map<pair<int, int>, int> m2;
            bedcounts1.insert(make_pair(*bc, m2));
            bedcounts2.insert(make_pair(*bc, m2));
        }
    }
    
    // Create scope to free bam reader
    {
        // Init BAM reader
        bam_reader reader(bamfile);
        if (has_out1){
            reader.write_header(out1);
        }        
        if (has_out2){
            reader.write_header(out2);
        }
        if (cell_barcode){
            reader.set_cb();
        }
        set<int> ahm_allkeys;
        if (alleleHMM){
            if (cell_barcode){
                fprintf(stdout, "chrm\tsnppos\tallelecount1\tallelecount2\tcell_barcode\n");
            }
            else{
                fprintf(stdout, "chrm\tsnppos\tallelecount1\tallelecount2\n");
        
            }
        }

        map<int, char>::iterator nextsnp1;
        map<int, char>::iterator nextsnp2;
        string curchrom = "";
        int winstart = 0;

        while(reader.next()){
            // Only count reads that are unmapped. If considering cell barcodes,
            // also only count reads with valid/passing cell barcodes.
            if (!reader.unmapped() && (!cell_barcode || 
                (reader.has_cb_z && all_barcodes.find(reader.cb_z) != all_barcodes.end()))){
                ++tot;

                char* chrom = reader.ref_id();
                if (strcmp(chrom, curchrom.c_str()) != 0){
                    if (alleleHMM && ahm_allkeys.size() > 0){
                        // Print and delete everything in it.
                        print_alleleHMM(curchrom, counts1, counts2, ahm_allkeys, -1, cell_barcode);
                    }
                    curchrom = chrom;
                    nextsnp1 = snps1[curchrom].begin();
                    nextsnp2 = snps2[curchrom].begin();
                    winstart = 0;
                }
                else if (alleleHMM){
                    // Print and delete everything in alleleHMM data structures before the given 
                    // start map coordinate
                    print_alleleHMM(curchrom, counts1, counts2, ahm_allkeys, reader.reference_start+1, 
                        cell_barcode);
                }
                // Catch up SNPs to the current position.
                while (nextsnp1 != snps1[curchrom].end() && nextsnp1->first < reader.reference_start + 1){
                    ++nextsnp1;
                }
                while (nextsnp2 != snps2[curchrom].end() && nextsnp2->first < reader.reference_start + 1){
                    ++nextsnp2;
                }
                // Check all SNPs that might fall within this read.
                map<int, char>::iterator read_snp1 = nextsnp1;
                map<int, char>::iterator read_snp2 = nextsnp2;
                // Count how many alleles the read contains that match version 1 and version 2
                // of each SNP
                int count_match1 = 0;
                int count_match2 = 0;
                while (read_snp1 != snps1[curchrom].end() && read_snp1->first >= reader.reference_start + 1 &&
                    read_snp1->first <= reader.reference_end){
                    char allele = reader.get_base_at(read_snp1->first);
                    //fprintf(stderr, "%c | 1 %c\n", allele, read_snp1->second);
                    //fprintf(stderr, "%s\t%d\t%c\n", curchrom.c_str(), read_snp1->first, allele);
                    if (allele == read_snp1->second){
                        count_match1++;
                        if (alleleHMM){
                            string bc = "";
                            if (cell_barcode && reader.has_cb_z){
                                bc = reader.cb_z;
                            }
                            if (counts1[bc].count(read_snp1->first) == 0){
                                counts1[bc].insert(make_pair(read_snp1->first, 1));
                            }
                            else{
                                counts1[bc][read_snp1->first]++;
                            }
                            ahm_allkeys.insert(read_snp1->first);
                        }
                    }
                    ++read_snp1;
                }
                while (read_snp2 != snps2[curchrom].end() && read_snp2->first >= reader.reference_start + 1 &&
                    read_snp2->first <= reader.reference_end){
                    char allele = reader.get_base_at(read_snp2->first);
                    //fprintf(stderr, "%c | 2 %c\n", allele, read_snp2->second);
                    //fprintf(stderr, "%s\t%d\t%c\n", curchrom.c_str(), read_snp2->first, allele);
                    if (allele == read_snp2->second){
                        count_match2++;
                        if (alleleHMM){
                            string bc = "";
                            if (cell_barcode && reader.has_cb_z){
                                bc = reader.cb_z;
                            }
                            if (counts2[bc].count(read_snp2->first) == 0){
                                counts2[bc].insert(make_pair(read_snp2->first, 1));
                            }
                            else{
                                counts2[bc][read_snp2->first]++;
                            }
                            ahm_allkeys.insert(read_snp2->first);
                        }
                    }
                    ++read_snp2;
                }
                if (paired){
                    // Store this for later.
                    char* ridbuf = reader.read_id();
                    string rid = ridbuf;
                    if (count_match1 > 0){
                        if (read_count_snp1.count(rid) == 0){
                            read_count_snp1.insert(make_pair(rid, count_match1));
                        }
                        else{
                            read_count_snp1[rid] += count_match1;
                        }
                    }
                    if (count_match2 > 0){
                        if (read_count_snp2.count(rid) == 0){
                            read_count_snp2.insert(make_pair(rid, count_match2));
                        }
                        else{
                            read_count_snp2[rid] += count_match2;
                        }
                    }
                }
                else{
                    //fprintf(stderr, "%d %d\n", count_match1, count_match2);
                    if (count_match1 > 0 && count_match2 == 0){
                        count_1++;
                        if (has_out1){
                            reader.write_record(out1);
                        }
                    }
                    else if (count_match2 > 0 && count_match1 == 0){
                        count_2++;
                        if (has_out2){
                            reader.write_record(out2);
                        }
                    }
                    else if (count_match1 > 0 && count_match2 > 0){
                        // ???
                       // fprintf(stderr, "%d %d\n", count_match1, count_match2);
                        count_chimeric++;
                    }
                    else{
                        // No informative SNPs
                        count_none++;
                    }
                }
            }
            if (tot % 10000 == 0){
                fprintf(stderr, "Processed %ld reads\r", tot);
            }
        }
        fprintf(stderr, "\n");
        if (alleleHMM && curchrom != ""){
            // Print everything counted on the final chromosome
            print_alleleHMM(curchrom, counts1, counts2, ahm_allkeys, -1, cell_barcode);
        }
    }
    if (paired && !alleleHMM){
        
        // Make a second pass.
        bam_reader reader(bamfile);
        if (cell_barcode){
            reader.set_cb();
        }
        fprintf(stderr, "Making second pass to catch mates of allele-specific reads\n");
        tot = 0;
        int max_winstart = 0;
        string curchrom = "";
        set<int> winstarts;
        
        set<pair<int, int> >::iterator curbed;
        bool bed_thischrom = false;

        while(reader.next()){
            if (!reader.unmapped()){
                ++tot;
                char* ridbuf = reader.read_id();
                string rid = ridbuf;
                char* chrom = reader.ref_id();
                if (strcmp(chrom, curchrom.c_str()) != 0){
                    if (window != -1){
                        // Process windows
                        max_winstart = 0;
                        curchrom = chrom;
                        process_windows(counts1, counts2, all_barcodes, curchrom, -1, 
                            window, offset, winstarts, cell_barcode);
                        winstarts.insert(0);
                    }
                    if (bedfile_given){
                        if (bed.count(chrom) > 0){
                            curbed = bed[chrom].begin();
                            bed_thischrom = true;
                        }
                        else{
                            bed_thischrom = false;
                        }                        
                    }
                }
                if (window != -1){
                    while (reader.reference_start > max_winstart + offset){
                        max_winstart += offset;
                        winstarts.insert(max_winstart);
                    }
                }
                if (window != -1 && *winstarts.begin() + window < reader.reference_start){
                    // Process windows out of range.
                    process_windows(counts1, counts2, all_barcodes, curchrom, reader.reference_start, window, 
                        offset, winstarts, cell_barcode);
                    
                }
                if (bedfile_given && bed_thischrom){
                    process_bed(bedcounts1, bedcounts2, all_barcodes, curchrom, reader.reference_start+1,
                        bed, cell_barcode);
                }
                if (read_count_snp1.count(rid) > 0){
                    if (read_count_snp2.count(rid) > 0){
                        count_chimeric++;
                    }
                    else{
                        count_1++;
                        if (window != -1){
                            string bc = "";
                            if (cell_barcode && reader.has_cb_z){
                                bc = reader.cb_z;
                            }
                            for (set<int>::iterator winstart = winstarts.begin(); 
                                winstart != winstarts.end(); ++winstart){
                                if (counts1[bc].count(*winstart) == 0){
                                    counts1[bc].insert(make_pair(*winstart, 1));
                                }
                                else{
                                    counts1[bc][*winstart]++;
                                }
                            }
                        }
                        if (bedfile_given && bed_thischrom && (!cell_barcode || reader.has_cb_z)){
                            set<pair<int, int> >::iterator bed_it2 = curbed;
                            while (bed_it2 != bed[curchrom].end() ){
                                if (bed_it2->first > reader.reference_start+1){
                                    break;
                                }   
                                else if (reader.reference_start+1 >= bed_it2->first && 
                                    reader.reference_start+1 <= bed_it2->second){
                                    string bc = "";
                                    if (cell_barcode){
                                        bc = reader.cb_z;
                                    }
                                    if (bedcounts1.count(bc) == 0){
                                        map<pair<int, int>, int> m;
                                        bedcounts1.insert(make_pair(bc, m));
                                    }
                                    if (bedcounts1[bc].count(*bed_it2) == 0){
                                        bedcounts1[bc].insert(make_pair(*bed_it2, 1));
                                    }
                                    else{
                                        bedcounts1[bc][*bed_it2]++;
                                    }
                                }
                                ++bed_it2;
                            }
                        }
                        if (has_out1){
                            reader.write_record(out1);
                        }
                    }
                }
                else if (read_count_snp2.count(rid) > 0){
                    count_2++;
                    if (window != -1){
                        string bc = "";
                        if (cell_barcode && reader.has_cb_z){
                            bc = reader.cb_z;
                        }
                        for (set<int>::iterator winstart = winstarts.begin();
                            winstart != winstarts.end(); ++winstart){
                            if (counts2[bc].count(*winstart) == 0){
                                counts2[bc].insert(make_pair(*winstart, 1));
                            }
                            else{
                                counts2[bc][*winstart]++;
                            }
                        }
                    }
                    if (bedfile_given && bed_thischrom && (!cell_barcode || reader.has_cb_z)){
                        set<pair<int, int> >::iterator bed_it2 = curbed;
                        while (bed_it2 != bed[curchrom].end() ){
                            if (bed_it2->first > reader.reference_start+1){
                                break;
                            }   
                            else if (reader.reference_start+1 >= bed_it2->first && 
                                reader.reference_start+1 <= bed_it2->second){
                                string bc = "";
                                if (cell_barcode){
                                    bc = reader.cb_z;
                                }
                                if (bedcounts2.count(bc) == 0){
                                    map<pair<int, int>, int> m;
                                    bedcounts2.insert(make_pair(bc, m));
                                }
                                if (bedcounts2[bc].count(*bed_it2) == 0){
                                    bedcounts2[bc].insert(make_pair(*bed_it2, 1));
                                }
                                else{
                                    bedcounts2[bc][*bed_it2]++;
                                }
                            }
                            ++bed_it2;
                        }
                    }

                    if (has_out2){
                        reader.write_record(out2);
                    }
                }
                else{
                    count_none++;
                }
            }   
            if (tot % 10000 == 0){
                fprintf(stderr, "Processed %ld reads\r", tot);
            }    
        }
        fprintf(stderr, "\n");
        if (window != -1){
            process_windows(counts1, counts2, all_barcodes, curchrom, -1, window, offset, winstarts, cell_barcode);  
        }
    }
    
    if (has_out1){
        bgzf_close(out1);
    }
    if (has_out2){
        bgzf_close(out2);
    }
    fprintf(stderr, "%ld total reads\n", tot);
    fprintf(stderr, "%ld matched allele 1 (%.2f)%%\n", count_1, 100*(float)count_1/(float)tot);
    fprintf(stderr, "%ld matched allele 2 (%.2f)%%\n", count_2, 100*(float)count_2/(float)tot);
    fprintf(stderr, "%ld chimeric (%.2f)%%\n", count_chimeric, 100*(float)count_chimeric/(float)tot);
    fprintf(stderr, "%ld lacked informative SNPs (%.2f)%%\n", count_none, 100*(float)count_none/(float)tot);
    return 0;
}

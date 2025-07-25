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
#include <cstdlib>
#include <utility>
#include <htslib/kseq.h>

using std::cout;
using std::endl;
using namespace std;

KSEQ_INIT(gzFile, gzread);
 
/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "fagrep fileIdx searchstr outdir file1.f(a|q)(.gz) file2.f(a|q)(.gz) file3.fa(a|q)(.gz) ...\n");
    fprintf(stderr, "Given one or more FAST(A/Q) files (paired, if multiple -- i.e. R1, R2, I1, I2), \
and a search sequence, writes out all reads containing the search sequence to output files, \
preserving pair info.\n");
    fprintf(stderr, "Positional arguments:\n");
    fprintf(stderr, "   fileIdx 1-based index of file to grep in (separate multiple with commas, no spaces)\n");
    fprintf(stderr, "   searchstr the sequence to grep for. If the sequence should come after n fixed bases at \
the beginning of a sequence, precede it by that many Ns. If it should come before n fixed bases at the end of a \
sequence, follow it by that many Ns. In other words, if searching for an 8-base UMI at the beginning of a read \
followed by \"ACGTAGT\", search for NNNNNNNNACGTAGT.\n");
fprintf(stderr, "If searching multiple files (fileIdx is comma separated), you must provide one search string per \
file index, also comma-separated, in the same order.\n");
    fprintf(stderr, "   outdir a directory to store output files, which will have the same name as input files\n");
    fprintf(stderr, "   file1 file2... one or more fasta or fastq files, optionally gzipped\n");
    exit(code);
}

void write_seqs(vector<kseq_t*>& seqs_in,
    vector<gzFile>& fp_out){
    
    for (int i = 0 ; i < seqs_in.size(); ++i){
        if (seqs_in[i]->qual.l > 0){
            // FASTQ
            gzwrite(fp_out[i], "@", 1);
            gzwrite(fp_out[i], seqs_in[i]->name.s, seqs_in[i]->name.l);
            gzwrite(fp_out[i], "\n", 1);
            gzwrite(fp_out[i], seqs_in[i]->seq.s, seqs_in[i]->seq.l);
            gzwrite(fp_out[i], "\n+\n", 3);
            gzwrite(fp_out[i], seqs_in[i]->qual.s, seqs_in[i]->qual.l);
            gzwrite(fp_out[i], "\n", 1);
        }
        else{
            // FASTA
            gzwrite(fp_out[i], ">", 1);
            gzwrite(fp_out[i], seqs_in[i]->name.s, seqs_in[i]->name.l);
            gzwrite(fp_out[i], "\n", 1);
            gzwrite(fp_out[i], seqs_in[i]->seq.s, seqs_in[i]->seq.l);
            gzwrite(fp_out[i], "\n", 1);
        }
    }
}

bool file_exists (const string& name) {
    ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char *argv[]) {    
    
    if (argc < 5){
        help(1);
    }    
    
    vector<int> file_indices;
    string file_indices_str = argv[1];
    size_t comma_pos = file_indices_str.find(",");
    while(comma_pos != string::npos){
        string field = file_indices_str.substr(0, comma_pos);
        file_indices_str = file_indices_str.substr(comma_pos + 1, file_indices_str.length()-comma_pos-1);
        int this_fileIdx = atoi(field.c_str());
        file_indices.push_back(this_fileIdx);
        comma_pos = file_indices_str.find(",");
    }
    file_indices.push_back(atoi(file_indices_str.c_str()));
    for (int x = 0; x < file_indices.size(); ++x){
        if (file_indices[x] <= 0 || file_indices[x] > argc-4){
            fprintf(stderr, "ERROR: invalid file index %d for %d input files\n", file_indices[x], argc-4);
            exit(1);
        }
        // Convert to 0-based.
        file_indices[x]--;
    }
    vector<string> searchstrs;
    string searchstr_all = argv[2];
    comma_pos = searchstr_all.find(",");
    while (comma_pos != string::npos){
        string field = searchstr_all.substr(0, comma_pos);
        searchstr_all = searchstr_all.substr(comma_pos + 1, searchstr_all.length()-comma_pos-1);
        
        searchstrs.push_back(field);
        comma_pos = searchstr_all.find(",");
    }
    searchstrs.push_back(searchstr_all);
    for (int x = 0; x < searchstrs.size(); ++x){
        for (int i = 0; i < searchstrs[x].length(); ++i){
            if (searchstrs[x][i] != 'A' && 
                searchstrs[x][i] != 'C' && 
                searchstrs[x][i] != 'G' && 
                searchstrs[x][i] != 'T' && 
                searchstrs[x][i] != 'N'){
                fprintf(stderr, 
                    "ERROR: illegal character %c detected in search string %s\n", 
                    searchstrs[x][i], searchstrs[x].c_str());
                exit(1);
            }
        }
    }
    if (searchstrs.size() != file_indices.size()){
        fprintf(stderr, "ERROR: %ld file indices and %ld search strings. Numbers of both must match.\n",
            file_indices.size(), searchstrs.size());
        exit(1);
    }
    
    vector<int> fixed_begins;
    vector<int> fixed_ends;
    for (int s_idx = 0; s_idx < searchstrs.size(); ++s_idx){
        int fixed_begin = -1;
        int fixed_end = -1;
        string searchstr = searchstrs[s_idx];
        for (int i = 0; i < searchstr.length(); ++i){
            if (searchstr[i] != 'N'){
                break;
            }
            else{
                fixed_begin = i + 1;
            }
        }
        for (int i = searchstr.length()-1; i >= 0; --i){
            if (searchstr[i] != 'N'){
                break;
            }
            else{
                if (fixed_end == -1){
                    fixed_end = 0;
                }
                fixed_end++;
            }
        }
        if (fixed_begin != -1 && fixed_end != -1){
            fprintf(stderr, "ERROR: runs of N are allowed at the beginning or end of a search string, but not both.\n");
            exit(1);
        }
        if (fixed_begin != -1){
            searchstr = searchstr.substr(fixed_begin, searchstr.length()-fixed_begin);
            searchstrs[s_idx] = searchstr;
        }
        if (fixed_end != -1){
            searchstr = searchstr.substr(0, searchstr.length()-fixed_end);
            searchstrs[s_idx] = searchstr;
        }
        fixed_begins.push_back(fixed_begin);
        fixed_ends.push_back(fixed_end);
    } 
    
    string outdir = argv[3];
    if (outdir == "." || outdir == "./"){
        outdir = "";
    }
    else if (outdir[outdir.length()-1] != '/'){
        outdir += "/";
    }

    vector<gzFile> fp_in;    
    vector<kseq_t*> seq_in;
    vector<gzFile> fp_out;
    vector<int> progress;

    // Parse input && ensure new file names don't coincide with input files
    bool approval_given = false;

    for (int i = 4; i < argc; ++i){
        string inpath = argv[i];
        string in_base = inpath.substr(inpath.find_last_of("/\\") + 1);
        string outpath = outdir + in_base;
        if (!approval_given && file_exists(outpath)){
            cerr << "Existing files found in output directory. Overwrite? (Y/N) ";
            char response;
            cin >> response;
            while (response != 'Y' && response != 'y' && response != 'N' &&
                response != 'n'){
                cerr << "Invalid response: " << response;
                cerr << "\nExisting files found in output directory. Overwrite? (Y/N) ";
                cin >> response;
            }
            if (response == 'n' || response == 'N'){
                // bail out.
                return 0;
            }
            else if (response == 'y' || response == 'Y'){
                approval_given = true;
            }
        }        
        
        gzFile fp_in_this = gzopen(inpath.c_str(), "r");
        gzFile fp_out_this = gzopen(outpath.c_str(), "w");
        
        fp_in.push_back(fp_in_this);
        fp_out.push_back(fp_out_this);

        kseq_t* seq_in_this = kseq_init(fp_in_this);
        
        seq_in.push_back(seq_in_this);
        
        progress.push_back(0);
    }
    
    long int matches = 0;
    long int reads = 0;

    while ((progress[0] = kseq_read(seq_in[0])) >= 0){
        
        // Advance everybody else
        for (int i = 1 ; i < seq_in.size(); ++i){
            progress[i] = kseq_read(seq_in[i]);
            if (progress[i] < 0 && progress[0] >= 0){
                fprintf(stderr, "ERROR: mismatched number of sequences\n");
                exit(1);
            }
            if (strcmp(seq_in[i]->name.s, seq_in[0]->name.s) != 0){
                fprintf(stderr, "ERROR: name mismatch in input sequences: %s and %s\n", seq_in[i]->name.s, seq_in[0]->name.s);
                exit(1);
            }
        }
        
        ++reads;

        // Need to find all sequences in all files to pass.
        int foundcount = 0;

        for (int grep_idx = 0; grep_idx != searchstrs.size(); ++grep_idx){
            string searchstr = searchstrs[grep_idx];
            int fileIdx = file_indices[grep_idx];
            int fixed_begin = fixed_begins[grep_idx];
            int fixed_end = fixed_ends[grep_idx];
            
            // Now grep appropriate sequence.
            string seq = seq_in[fileIdx]->seq.s;
    
            int searchIdx = seq.find(searchstr);
            if (searchIdx != string::npos){
                if ((fixed_begin == -1 || searchIdx == fixed_begin) &&
                    (fixed_end == -1 || searchIdx == seq.length()-fixed_end-searchstr.length())){
                    ++foundcount;
                }
            }
        }

        if (foundcount == searchstrs.size()){
            // Found all necessary sequences in all files.
            ++matches;
            if (matches % 1000 == 0){
                fprintf(stderr, "%ld matches in %ld reads (%.1f%%)\r", 
                    matches, reads, 100.0 * ((float)matches/(float)reads));
            }
            write_seqs(seq_in, fp_out);
        }
    }
    fprintf(stderr, "\n");
    
    for (int i = 0; i < seq_in.size(); ++i){
        kseq_destroy(seq_in[i]);
        gzclose(fp_out[i]);
    }
    
    return 0;
}

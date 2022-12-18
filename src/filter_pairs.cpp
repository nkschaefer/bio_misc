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
    fprintf(stderr, "filter_pairs [OPTIONS]\n");
    fprintf(stderr, "Given two FASTA/FASTQ files (IN READ-ID SORTED ORDER), creates two \
new FASTA/FASTQ files containing only paired reads in the appropriate order.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --forward -1 The forward read file.\n");
    fprintf(stderr, "    --reverse -2 The reverse read file.\n");
    fprintf(stderr, "    --out1 -3 The output file for forward reads");
    fprintf(stderr, "    --out2 -4 The output file for reverse reads\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

void write_fastq(kseq_t* seq, gzFile& out){
    int buflen = seq->name.l + 3;
    if (seq->comment.l > 0){
        buflen += seq->comment.l + 1;
    }
    char buf[buflen];
    if (seq->comment.l > 0){
        sprintf(buf, "@%s %s\n", seq->name.s, seq->comment.s);
    }
    else{
        sprintf(buf, "@%s\n", seq->name.s);
    }
    gzwrite(out, buf, buflen-1);
    char buf2[seq->seq.l + 1];
    sprintf(buf2, "%s\n", seq->seq.s);
    gzwrite(out, buf2, seq->seq.l + 1);
    char buf3[3];
    sprintf(buf3, "+\n");
    gzwrite(out, buf3, 2);
    sprintf(buf2, "%s\n", seq->qual.s);
    gzwrite(out, buf2, seq->qual.l + 1);
}

void write_fasta(kseq_t* seq, gzFile& out){
    int buflen = seq->name.l + 3;
    if (seq->comment.l > 0){
        buflen += seq->comment.l + 1;
    }
    char buf[buflen];
    if (seq->comment.l > 0){
        sprintf(buf, ">%s %s\n", seq->name.s, seq->comment.s);
    }
    else{
        sprintf(buf, ">%s\n", seq->name.s);
    }
    gzwrite(out, buf, buflen-1);
    char buf2[seq->seq.l + 1];
    sprintf(buf2, "%s\n", seq->seq.s);
    gzwrite(out, buf2, seq->seq.l + 1);
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
       {"forward", required_argument, 0, '1'},
       {"reverse", required_argument, 0, '2'},
       {"out1", required_argument, 0, '3'},
       {"out2", required_argument, 0, '4'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string forwardfile;
    string revfile;
    string out1file;
    string out2file;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "1:2:3:4:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case '1':
                forwardfile = optarg;
                break;
            case '2':
                revfile = optarg;
                break;
            case '3':
                out1file = optarg;
                break;
            case '4':
                out2file = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    if (forwardfile.length() == 0 || revfile.length() == 0 || out1file.length() == 0 ||
        out2file.length() == 0){
        fprintf(stderr, "ERROR: one or more required arguments was missing.\n");
        exit(1);
    }
    
    int f_progress;
    int r_progress;
    
    gzFile f_fp = gzopen(forwardfile.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening file %s for reading.\n", forwardfile.c_str());
        exit(1);
    }
    gzFile r_fp = gzopen(revfile.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening file %s for reading.\n", revfile.c_str());
        exit(1);
    }
    
    if (out1file.substr(out1file.length()-3, 3) != ".gz"){
        out1file += ".gz";
    }
    if (out2file.substr(out2file.length()-3, 3) != ".gz"){
        out2file += ".gz";
    }
    
    gzFile f_out_fp = gzopen(out1file.c_str(), "w");
    if (!f_out_fp){
        fprintf(stderr, "ERROR opening file %s for writing.\n", out1file.c_str());
        exit(1);
    }
    gzFile r_out_fp = gzopen(out2file.c_str(), "w");
    if (!r_out_fp){
        fprintf(stderr, "ERROR opening file %s for writing.\n", out2file.c_str());
        exit(1);
    }
    
    kseq_t* seq_f = kseq_init(f_fp);
    kseq_t* seq_r = kseq_init(r_fp);
    
    long int seq_count = 0;
    long int progress = 1000;
    long int count_unpaired = 0;
    
    // Initial read
    f_progress = kseq_read(seq_f);
    r_progress = kseq_read(seq_r);
        
    while (f_progress >= 0 && r_progress >= 0){
        // If IDs don't match, fast forward the lagging file until we get a match or
        // go past the other one.
        string id_f = seq_f->name.s;
        string id_r = seq_r->name.s;
        while (id_f != id_r){
            while (id_f < id_r){
                f_progress = kseq_read(seq_f);
                count_unpaired++;
                if (f_progress < 0){
                    break;
                }
                id_f = seq_f->name.s;
            }
            if (f_progress < 0){
                break;
            }
            while (id_r < id_f){
                r_progress = kseq_read(seq_r);
                count_unpaired++;
                if (r_progress < 0){
                    break;
                }
                id_r = seq_r->name.s;
            }
            if (r_progress < 0){
                break;
            }
        }
        if (f_progress >= 0 && r_progress >= 0 && id_f == id_r){
            // Print sequences.
            if (seq_f->qual.l > 0){
                // FASTQ
                write_fastq(seq_f, f_out_fp);
            }
            else{
                write_fasta(seq_f, f_out_fp);
            }
            if (seq_r->qual.l > 0){
                // FASTQ
                write_fastq(seq_r, r_out_fp);
            }
            else{
                write_fasta(seq_r, r_out_fp);
            }
            seq_count++;
            if (seq_count % progress == 0){
                fprintf(stderr, "Printed %ld pairs\r", seq_count);
            }
            f_progress = kseq_read(seq_f);
            r_progress = kseq_read(seq_r);
        }
    }
    
    fprintf(stderr, "Printed %ld pairs\n", seq_count);
    fprintf(stderr, "Removed %ld unpaired reads\n", count_unpaired);
    
    // Clean up.
    gzclose(f_out_fp);
    gzclose(r_out_fp);
    
    kseq_destroy(seq_f);
    kseq_destroy(seq_r);
    
    return 0;
}

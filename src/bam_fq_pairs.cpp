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
#define MAX_ID_LEN (101)
#define MAX_SEQ_LEN (501)
#define MAX_BC_LEN (51)

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bam_fq_pairs [OPTIONS]\n");
    fprintf(stderr, "Replacement for samtools fastq. Expects BAM input sorted by read name. \
Pairs reads (and outputs unpaired reads separately). Also dumps all tags from each BAM entry \
into the FASTQ comment field. These can then be re-inserted into a BAM using bwa mem -C.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --r1 -1 The output file for forward, paired reads\n");
    fprintf(stderr, "    --r2 -2 The output file for reverse, paired reads\n");
    fprintf(stderr, "    --single -s The output file for unpaired reads\n");
    fprintf(stderr, "    --scRNA -r Specify if the BAM file is 10x scRNA-seq data. In that case, \
this program will create forward reads consisting of cell barcode sequences.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

/**
 * Print FASTQ data to the output file.
 */
void write_fq(gzFile& out, const char* id, const char* comment, const char* seq, const char* qual){
    gzwrite(out, "@", 1);
    gzwrite(out, id, strlen(id));
    gzwrite(out, " ", 1);
    gzwrite(out, comment, strlen(comment));
    gzwrite(out, "\n", 1);
    gzwrite(out, seq, strlen(seq));
    gzwrite(out, "\n+\n", 3);
    gzwrite(out, qual, strlen(qual));
    gzwrite(out, "\n", 1);
}

/*
// Store all 10X barcodes.
struct bcs_10x{
    char bc[MAX_BC_LEN];
    char bx[MAX_BC_LEN];
    char st[MAX_BC_LEN];
    char rx[MAX_BC_LEN];
    char qx[MAX_BC_LEN];
    char tr[MAX_BC_LEN];
    char tq[MAX_BC_LEN];
    inline void populate(bam_reader& reader){
        if (reader.has_bc_z){
            sprintf(&this->bc[0], "%s", reader.bc_z);
        }
        else{
            this->bc[0] = '\0';
        }
        if (reader.has_bx_z){
            sprintf(&this->bx[0], "%s", reader.bx_z);
        }
        else{
            this->bx[0] = '\0';
        }
        if (reader.has_st_z){
            sprintf(&this->st[0], "%s", reader.st_z);
        }
        else{
            this->st[0] = '\0';
        }
        if (reader.has_rx_z){
            sprintf(&this->rx[0], "%s", reader.rx_z);
        }
        else{
            this->rx[0] = '\0';
        }
        if (reader.has_qx_z){
            sprintf(&this->qx[0], "%s", reader.qx_z);
        }
        else{
            this->qx[0] = '\0';
        }
        if (reader.has_tr_z){
            sprintf(&this->tr[0], "%s", reader.tr_z);
        }
        else{
            this->tr[0] = '\0';
        }
        if (reader.has_tq_z){
            sprintf(&this->tq[0], "%s", reader.tq_z);
        }
        else{
            this->tq[0] = '\0';
        }
    };
};
*/

void parse_cur_read(bam_reader& reader, char* id, char* comment, 
    char* seq, char* qual, char* rtype){
    
    char* cur_id = reader.read_id();
    sprintf(id, "%s", cur_id);
    reader.get_seq(seq);
    reader.get_qual(qual);
    if (reader.read1()){
        rtype[0] = '1';
    }
    else{
        rtype[0] = '2';
    }
    
    int comment_pos = 0;
    static vector<string> tags = {"CR", "CY", "CB", "RG", "BX"};
    for (int i = 0; i < tags.size(); ++i){
        uint8_t* aux_bin = bam_aux_get(reader.reader, tags[i].c_str());
        if (aux_bin != NULL){
            char* val = bam_aux2Z(aux_bin);
            if (comment_pos > 0){
                sprintf(comment + comment_pos, "\t%s:Z:%s", tags[i].c_str(), val);
                comment_pos++;
            }
            else{
                sprintf(comment + comment_pos, "%s:Z:%s", tags[i].c_str(), val);
            }
            comment_pos += 5 + strlen(val);
        } 
    }
    /*
    // Load aux data
    uint8_t* auxdat = bam_get_aux(reader.reader);
    vector<string> auxtags;
    char tagname[2];
    
    while (true){
        tagname[0] = auxdat[0];
        tagname[1] = auxdat[1];
        switch(auxdat[3]){
            case 
        }
        auxdat += 3;
        // Read data

    } 
    fprintf(stderr, "%c%c%c%c%c\n", auxdat[0], auxdat[1], auxdat[2], auxdat[3], auxdat[4]);
    */
    /*
    sprintf(comment, "%s:N:0:BC:Z:%s:ST:Z:%s_RX:Z:%s_QX:Z:%s_TR:Z:%s_TQ:Z:%s", \
        rtype, &bcs.bc[0], &bcs.st[0], &bcs.rx[0], &bcs.qx[0], &bcs.tr[0], &bcs.tq[0]);
    */
    // Per https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines:
    // barcodes embedded as a coment of the form: BX:Z:ACTCGACTGACTAGCT-1
    //sprintf(comment, "BX:Z:%s", &bcs.bx[0]);
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
       {"r1", required_argument, 0, '1'},
       {"r2", required_argument, 0, '2'},
       {"single", required_argument, 0, 's'},
       {"scRNA", no_argument, 0, 'r'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string r1file;
    string r2file;
    bool has_r1 = false;
    bool has_r2 = false;
    bool has_paired = false;
    string sfile;
    bool has_single = false;
    bool scRNA = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:1:2:s:rh", long_options, &option_index )) != -1){
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
                r1file = optarg;
                has_r1 = true;
                break;
            case '2':
                r2file = optarg;
                has_r2 = true;
                break;
            case 's':
                sfile = optarg;
                has_single = true;
                break;
            case 'r':
                scRNA = true;
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
    has_paired = has_r1 && has_r2;
    if (!has_paired && !has_single){
        fprintf(stderr, "ERROR: at least one of either -1 and -2 or -s must be provided.\n");
        exit(1);
    }
    
    gzFile r1;
    gzFile r2;
    if (has_paired){
        if (r1file.length() < 3 || r1file.substr(r1file.length()-3, 3) != ".gz"){
            r1file += ".gz";
        }
        if (r2file.length() < 3 || r2file.substr(r2file.length()-3, 3) != ".gz"){
            r2file += ".gz";
        }
        r1 = gzopen(r1file.c_str(), "w");
        if (r1 == Z_NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", r1file.c_str());
            exit(1);
        }
        r2 = gzopen(r2file.c_str(), "w");
        if (r2 == Z_NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", r2file.c_str());
            exit(1);
        }
    }
    gzFile s;
    if (has_single){
        if (sfile.length() < 3 || sfile.substr(sfile.length()-3, 3) != ".gz"){
            sfile += ".gz";
        }
        s = gzopen(sfile.c_str(), "w");
        if (s == Z_NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", sfile.c_str());
            exit(1);
        }
    }

    // Init BAM reader
    bam_reader reader(bamfile);
    reader.set_10x();
    if (scRNA){
        string bcstr = "CB";
        reader.set_bc_tag(bcstr);
    }

    char id[MAX_ID_LEN];
    char comment[5000];
    char seq[MAX_SEQ_LEN];
    char qual[MAX_SEQ_LEN];
    
    // If we need to write out cell barcodes as reads, fill in a fake
    // quality string the same length as barcode sequences
    char bc_qual[17];
    for (int i = 0; i < 16; ++i){
        bc_qual[i] = '?';
    } 
    bc_qual[16] = '\0';
    char bcbuf[100];

    bool has_prev = false;
    // Need to track this for 10X barcode formatting in FASTQ
    char rtype[2];
    rtype[1] = '\0';
    
    while(reader.next()){
        
        // If not tracking read pairs, just print everything right away.
        if (!has_paired){
            // Parse current sequence.
            parse_cur_read(reader, &id[0], &comment[0], &seq[0], &qual[0], &rtype[0]);
            write_fq(s, &id[0], &comment[0], &seq[0], &qual[0]);
        }
        else{
            if (reader.paired() && reader.read1()){
                if (has_prev){
                    // Must be unpaired.
                    if (has_single){
                        write_fq(s, &id[0], &comment[0], &seq[0], &qual[0]);
                    }
                }
                // Parse current sequence.
                parse_cur_read(reader, &id[0], &comment[0], &seq[0], &qual[0], &rtype[0]);
                has_prev = true;
            }
            else if (reader.paired() && reader.read2()){
                // See whether IDs match.
                char* cur_id = reader.read_id();
                bool had_r1 = false;
                if (strcmp(cur_id, &id[0]) == 0 && rtype[0] == '1'){
                    write_fq(r1, &id[0], &comment[0], &seq[0], &qual[0]);
                    had_r1 = true;
                }
                // Parse current sequence
                parse_cur_read(reader, &id[0], &comment[0], &seq[0], &qual[0], &rtype[0]);
                if (had_r1){
                    write_fq(r2, &id[0], &comment[0], &seq[0], &qual[0]);
                }
                else{
                    write_fq(s, &id[0], &comment[0], &seq[0], &qual[0]);
                }
                has_prev = false;
            }
            else{
                // Unpaired read. Just print it right away.
                parse_cur_read(reader, &id[0], &comment[0], &seq[0], &qual[0], &rtype[0]);
                if (scRNA && reader.bc.length() > 0){
                    bool terminated = false;
                    for (int i = 0; i < reader.bc.length(); ++i){
                        if (reader.bc[i] == 'A' || reader.bc[i] == 'C' || reader.bc[i] == 'G' || reader.bc[i] == 'T'){
                            bcbuf[i] = reader.bc[i];
                        }
                        else{
                            bcbuf[i] = '\0';
                            terminated = true;
                            break;
                        }
                    }    
                    if (!terminated){
                        bcbuf[reader.bc.length()] = '\0';
                    }               
                    write_fq(r1, &id[0], &comment[0], &bcbuf[0], &bc_qual[0]);
                    write_fq(r2, &id[0], &comment[0], &seq[0], &qual[0]); 
                }
                else{
                    write_fq(s, &id[0], &comment[0], &seq[0], &qual[0]);
                }
                has_prev = false;
            }
        }
    }
    
    if (has_paired){
        gzclose(r1);
        gzclose(r2);
    }
    if (has_single){
        gzclose(s);
    }
    
    return 0;
}

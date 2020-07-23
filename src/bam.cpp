/**
 * BAM reader wrapper
 */
#include <string> 
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "bam.h"

using namespace std;

// Define bases as represented in HTSLIB
uint8_t BASE_A = 1;
uint8_t BASE_C = 2;
uint8_t BASE_G = 4;
uint8_t BASE_T = 8;
uint8_t BASE_N = 15;

/**
 * Initialization function.
 *  @param bamfile
 *      A string representation of the name of the BAM file to parse
 */
bam_reader::bam_reader(string& bamfile){
    // Initialize BAM file reader
    this->filename = bamfile;
    this->fp = hts_open(filename.c_str(), "r");
    this->header = sam_hdr_read(fp);
    this->reader = bam_init1();
    this->has_bc_tag = false;
    this->prevtid = -1;
}

/**
 * Destructor
 */
bam_reader::~bam_reader(){
    bam_destroy1(this->reader);
    sam_close(this->fp);
}

/** 
 * Sets a barcode tag to look at.
 *  @param bctag
 *      A C++ string representation of the name of the tag of interest (i.e. BX)
 */
void bam_reader::set_bc_tag(string& bctag){
    this->bc_tag = bctag;
    this->has_bc_tag = true;
}

/**
 * Alternatively, set this to try to extract every 10X Genomics barcode in every read.
 */
void bam_reader::set_10x(){
    this->bcs_10x = true;
}

/**
 * Stop looking for every 10X Genomics barcode in every read.
 */
void bam_reader::unset_10x(){
    this->bcs_10x = false;
    
    this->has_bc_z = false;
    this->has_st_z = false;
    this->has_rx_z = false;
    this->has_qx_z = false;
    this->has_tr_z = false;
    this->has_tq_z = false;
}

/**
 * Get the next sequence from the BAM file.
 *
 * @return true if still reading (not EOF); false if finished (EOF).
 * 
 */
bool bam_reader::next(){
    while(sam_read1(this->fp, this->header, this->reader) > 0){
        if (this->reader->core.tid != this->prevtid){
            //fprintf(stderr, "Reading reference seq %s\n", this->header->target_name[this->reader->core.tid]);
            this->prevtid = this->reader->core.tid;
        }
        // Make insertion size positive
        this->isize = this->reader->core.isize;
        if (this->isize < 0){
            this->isize = 0 - this->isize;
        }
        if (this->has_bc_tag){
            // Retrieve barcode sequence.
            uint8_t* bc_bin = bam_aux_get(this->reader, this->bc_tag.c_str());
            if (bc_bin != NULL){
                // Convert to string
                char* bc_char = bam_aux2Z(bc_bin);
                // Convert to C++ style string
                this->bc = std::string(bc_char);
            }
        }
        if (this->bcs_10x){
            // Attempt to extract every 10X Genomics barcode tag.
            this->has_bc_z = false;
            this->has_bx_z = false;
            this->has_st_z = false;
            this->has_rx_z = false;
            this->has_qx_z = false;
            this->has_tr_z = false;
            this->has_tq_z = false;
            uint8_t* bc_bin = bam_aux_get(this->reader, "BC");
            if (bc_bin != NULL){
                this->has_bc_z = true;
                this->bc_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "BX");
            if (bc_bin != NULL){
                this->has_bx_z = true;
                this->bx_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "ST");
            if (bc_bin != NULL){
                this->has_st_z = true;
                this->st_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "RX");
            if (bc_bin != NULL){
                this->has_rx_z = true;
                this->rx_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "QX");
            if (bc_bin != NULL){
                this->has_qx_z = true;
                this->qx_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "TR");
            if (bc_bin != NULL){
                this->has_tr_z = true;
                this->tr_z = bam_aux2Z(bc_bin);
            }
            bc_bin = bam_aux_get(this->reader, "TQ");
            if (bc_bin != NULL){
                this->has_tq_z = true;
                this->tq_z = bam_aux2Z(bc_bin);
            }
        }
        return true;
    }
    return false; // EOF
}

// FLAG FUNCTIONS
// To avoid needing to constantly look these up, all are defined here.
// For reference, look at the function bam_flag2str() in sam.c in HTSLib.

bool bam_reader::paired(){
    return (this->reader->core.flag&BAM_FPAIRED);
}
bool bam_reader::proper_pair(){
    return (this->reader->core.flag&BAM_FPROPER_PAIR);
}
bool bam_reader::unmapped(){
    return (this->reader->core.flag&BAM_FUNMAP);
}
bool bam_reader::mate_unmapped(){
    return (this->reader->core.flag&BAM_FMUNMAP);
}
bool bam_reader::reverse(){
    return (this->reader->core.flag&BAM_FREVERSE);
}
bool bam_reader::mate_reverse(){
    return (this->reader->core.flag&BAM_FMREVERSE);
}
bool bam_reader::read1(){
    return (this->reader->core.flag&BAM_FREAD1);
}
bool bam_reader::read2(){
    return (this->reader->core.flag&BAM_FREAD2);
}
bool bam_reader::secondary(){
    return (this->reader->core.flag&BAM_FSECONDARY);
}
bool bam_reader::qcfail(){
    return (this->reader->core.flag&BAM_FQCFAIL);
}
bool bam_reader::dup(){
    return (this->reader->core.flag&BAM_FDUP);
}
bool bam_reader::supplementary(){
    return (this->reader->core.flag&BAM_FSUPPLEMENTARY);
}

/**
 * Returns the ID of the reference sequence for the current read.
 *  @return a C++ string representation of the reference genome sequence (chromosome) ID
 */
char* bam_reader::ref_id(){
    return this->header->target_name[this->reader->core.tid];
}

/**
 * Returns true if the read and mate are mapped to different chromosomes.
 */
bool bam_reader::chimeric(){
    if (! (this->reader->core.flag&BAM_FPAIRED)){
        // This isn't a paired read.
        return false;
    }
    return this->reader->core.tid != this->reader->core.mtid;
}

/**
 * Returns the ID of the current read
 *  @return a C++ string representation of the read ID
 */
char* bam_reader::read_id(){
    return bam_get_qname(this->reader);
}

/**
 * Retrieve the upstream-most (0-based ref genome coordinates) mapping position of the read
 *  @return 
 *      An integer coordinate of the first mapped position
 */
int bam_reader::map_start(){
    return this->reader->core.pos;
}

/**
 * Retrieve the downstream-most (0-based ref genome coordinates) mapping position of the read
 *  @return
 *      An integer coordinate of the last mapped position
 */
int bam_reader::map_end(){
    return bam_endpos(this->reader);
}

/**
 * Places a string representation of the read in provided char buffer.
 * It will be in reference sequence order/coordinates/complementarity.
 *
 * @param buf buffer in which to store the sequence.
 * @param int bases The number of bases to retrieve
 *      -1 or higher number than read length means retrieve entire sequence.
 */
void bam_reader::get_seq(char* seqbuf){
    // Fetch read sequence (binary format)
    uint8_t* seq_bin = bam_get_seq(this->reader);
    
    long int len = this->reader->core.l_qseq;

    seqbuf[len] = '\0';
    
    for (int i = 0; i < len; ++i){
        // Retrieve binary representation of base
        uint8_t base = bam_seqi(seq_bin, i);
        // Convert to character (if this is a reverse read, take its complement)
        char basechr;
        if (base == BASE_A){
            seqbuf[i] = 'A';
        }
        else if (base == BASE_C){
            seqbuf[i] = 'C';
        }
        else if (base == BASE_G){
            seqbuf[i] = 'G';
        }
        else if (base == BASE_T){
            seqbuf[i] = 'T';
        }
        else if (base == BASE_N){
            seqbuf[i] = 'N';
        }
        else{
            seqbuf[i] = 'N';
        }
    }
}

/**
 * Returns a char array pointer representation of the read's quality string.
 *
 */
void bam_reader::get_qual(char* buf){
    long int len = this->reader->core.l_qseq;
    uint8_t* buf_internal = bam_get_qual(this->reader);
    buf[len] = '\0';
    for (int i = 0; i < len; ++i){
        buf[i] = buf_internal[i] +33;
    }
    //sprintf(buf, "%s", buf_internal);
}

/**
 * Print the header to the given file in BAM format.
 *
 * @param fp the file to print to
 */
void bam_reader::write_header(BGZF* fp){
    int status = bam_hdr_write(fp, this->header);
}

/**
 * Print the current record to the given file in BAM format.
 */
void bam_reader::write_record(BGZF* fp){
    int status = bam_write1(fp, this->reader);
}
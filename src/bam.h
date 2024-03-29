#ifndef _BAMWRAPPER
#include <string>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

class bam_reader{
    private:
        std::string filename;
        samFile* fp;
        bam_hdr_t* header;
        std::string bc_tag;
        int prevtid;
        bool bcs_10x; // Do we want to try to extract every 10X Genomics barcode (see below)?
        bool cb; // Should we look for cell barcode tags for every entry?
        int32_t get_query_start();
        int32_t get_query_end();
        
    public:
        bam1_t* reader;
        std::string bc;
        
        long int reference_start;
        long int reference_end;
        long int reference_length;
        long int query_start;
        long int query_end;
        long int query_length;
        
        // Define stuff to store all potential 10X barcode data
        char* bc_z;
        char* bx_z;
        char* st_z;
        char* rx_z;
        char* qx_z;
        char* tr_z;
        char* tq_z;
        char* cb_z;

        bool has_bc_z;
        bool has_bx_z;
        bool has_st_z;
        bool has_rx_z;
        bool has_qx_z;
        bool has_tr_z;
        bool has_tq_z;
        bool has_cb_z;

        long int isize;
        bool has_bc_tag;
        
    // Constructor/destructor
    bam_reader(std::string&);
    ~bam_reader();
    
    void clear_read_groups_hdr();    
    void add_read_group_hdr(const std::string& id, 
        const std::string& sm, 
        const std::string& lib, 
        const std::string& pu, 
        const std::string& pl);
    void add_read_group_read(const std::string&);
    
    void add_string_tag(const std::string&, const std::string&);

    void set_bc_tag(std::string&);
    void set_cb();
    void set_10x();
    void unset_10x();
    
    bool next();
    
    // Flag functions
    bool paired();
    bool proper_pair();
    bool unmapped();
    bool mate_unmapped();
    bool reverse();
    bool mate_reverse();
    bool read1();
    bool read2();
    bool secondary();
    bool qcfail();
    bool dup();
    bool supplementary();
    
    char* ref_id();
    bool chimeric();
    char* read_id();
    int map_start();
    int map_end();
    void get_seq(char*);
    void get_qual(char*);
    void write_header(BGZF*);
    void write_record(BGZF*);
    char get_base_at(long int);
    
};

#define _BAMWRAPPER
#endif

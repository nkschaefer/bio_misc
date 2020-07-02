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
#include "treeNode.h"

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

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "eig_upgma [OPTIONS]\n");
    fprintf(stderr, "Given a genotype matrix (i.e. the .geno file from bcf2eigenstrat), \
computes pairwise distances between samples and creates a UPGMA tree.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --geno -g The input genotype matrix file.\n");
    fprintf(stderr, "    --samples -s The file listing sample names.\n");
    fprintf(stderr, "    --ploidy -p The file listing ploidy of each sample");
    fprintf(stderr, "    --distmat -d Write the distmat to this file\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

void parse_samples(string& filename, vector<string>& samplenames){
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
            samplenames.push_back(line);
        }
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

/**
 * Initialize a pre-declared distance matrix by setting all elements that will
 * be used and accessed to 0 and all other elements to -1.
 */
void init_distmat(vector<vector<float> >&dist_mat, int dim){
    for (int i = 0; i < dim-1; ++i){
        vector<float> row;
        for (int j = 0; j <= i; ++j){
            row.push_back(-1);
        }
        for (int j = i+1; j < dim; ++j){
            row.push_back(0);
        }
        dist_mat.push_back(row);
    }
}

void upgma(int num_haplotypes, 
    vector<vector<float> >& dist_mat,
    bool write_distmat, 
    string& distmatfile, 
    vector<string>& indvs){
    // Now we can use UPGMA to build a tree.
        
    // For each entry in the current distance matrix, store a pointer to a treeNode
    // that that distance matrix entry represents.
    
    fprintf(stderr, "Building tree\n");
    treeNode* tree = NULL;
    
    vector<treeNode*> dist_mat_nodes;
    vector<float> node_heights;
    for (int i = 0 ; i < num_haplotypes; ++i){
        cladeset leaves;
        leaves.set(num_haplotypes-i-1);
        treeNode* node = new treeNode(num_haplotypes, leaves);
        dist_mat_nodes.push_back(node);
        node_heights.push_back(0);
    }
    
    int dist_mat_dim = num_haplotypes;
    
    if (write_distmat){
        FILE* distmatf = fopen(distmatfile.c_str(), "w");
        if (indvs.size() == num_haplotypes){
            for (int i = 0; i < num_haplotypes; ++i){
                fprintf(distmatf, "\t%s", indvs[i].c_str());
            }
        }
        else{
            for (int i = 0; i < num_haplotypes; ++i){
                fprintf(distmatf, "\thap%d", i+1);
            }
        }
        fprintf(distmatf, "\n");
        for (int i = 0; i < num_haplotypes; ++i){            
            if (indvs.size() == num_haplotypes){
                fprintf(distmatf, "%s\t", indvs[i].c_str());
            }
            else{
                fprintf(distmatf, "hap%d\t", i+1);
            }
            // Lower values
            for (int j = 0; j < i; ++j){
                fprintf(distmatf, "\t%f", dist_mat[j][i]);
            }
            // Matching value
            fprintf(distmatf, "\t0");
            // Higher values
            for (int j = i+1; j < num_haplotypes; ++j){
                fprintf(distmatf, "\t%f", dist_mat[i][j]);
            }
            fprintf(distmatf, "\n");
        }
        fclose(distmatf);
    }
    
    while(dist_mat_dim > 1){
        fprintf(stderr, "Collapsing %d nodes...\r", dist_mat_dim);
        
        // Get the smallest distance in the matrix.
        float smallest_dist = dist_mat[0][1];
        int smallest_i = 0;
        int smallest_j = 1;
        for (int i = 0; i < dist_mat_dim-1; ++i){
            for (int j = i+1; j < dist_mat_dim; ++j){
                // Use -1 to represent things removed from the matrix.
                if (dist_mat[i][j] != -1 && dist_mat[i][j] < smallest_dist){
                    smallest_dist = dist_mat[i][j];
                    smallest_i = i;
                    smallest_j = j;
                }
            }
        }
        
        // Build a node connecting the closest things.
        treeNode* a = dist_mat_nodes[smallest_i];
        treeNode* b = dist_mat_nodes[smallest_j];
        cladeset a_subtree = a->subtree_leaves();
        cladeset b_subtree = b->subtree_leaves();
        cladeset parent_leaves;
        treeNode* ab = new treeNode(num_haplotypes, parent_leaves);
        a->parent = ab;
        b->parent = ab;
        
        a->dist = a->dist_norm = smallest_dist/2 + node_heights[smallest_j]/2;
        b->dist = b->dist_norm = smallest_dist/2 + node_heights[smallest_i]/2;
        
        ab->children.push_back(a);
        ab->children.push_back(b);
        
        int a_weight = a_subtree.count();
        int b_weight = b_subtree.count();
        
        float ab_height = (node_heights[smallest_i] + node_heights[smallest_j] + smallest_dist);
        
        // Rebuild distance matrix.
        dist_mat_dim--;
        
        // Create new distance matrix and matrix of tree pointers.
        // Allocate the same amount of space, so assignment continues to work,
        // but ignore the last elements.
        vector<vector<float> > dist_mat_new;
        init_distmat(dist_mat_new, dist_mat_dim);
        vector<treeNode*> dist_mat_nodes_new;
        vector<float> node_heights_new;
        
        // Replace previous i index with new node. Delete previous j index.
        for (int i = 0; i < dist_mat_dim; ++i){
            if (i != smallest_j){
                int new_i = i;
                if (i > smallest_j){
                    new_i--;
                }
                for (int j = i+1; j < dist_mat_dim+1; ++j){
                    if (j != smallest_j){
                        int new_j = j;
                        if (j > smallest_j){
                            new_j--;
                        }
                        if (i == smallest_i){
                            // Distance will be from the new node to the other stuff.
                            if (j > smallest_j){
                                dist_mat_new[new_i][new_j] = (dist_mat[smallest_i][j]*a_weight + dist_mat[smallest_j][j]*b_weight)/(a_weight+b_weight);
                            }
                            else{
                                dist_mat_new[new_i][new_j] = (dist_mat[smallest_i][j]*a_weight + dist_mat[j][smallest_j]*b_weight)/(a_weight+b_weight);
                            }
                        }
                        else{
                            // Copy over previous value.
                            dist_mat_new[new_i][new_j] = dist_mat[i][j];
                        }
                    }
                }
            }
        }
        
        // Update other stuff
        for (int i = 0; i < dist_mat_dim+1; ++i){
            if (i == smallest_i){
                // Add info for new node
                dist_mat_nodes_new.push_back(ab);
                node_heights_new.push_back(ab_height);
            }
            else if (i == smallest_j){
                // Skip index.
            }
            else{
                dist_mat_nodes_new.push_back(dist_mat_nodes[i]);
                node_heights_new.push_back(node_heights[i]);
            }
        }
        
        dist_mat = dist_mat_new;
        dist_mat_nodes = dist_mat_nodes_new;
        node_heights = node_heights_new;
        
    }

    tree = dist_mat_nodes[0];
    printf("%s\n", tree->newick(false, true, indvs).c_str());
    
    
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
       {"geno", required_argument, 0, 'g'},
       {"samples", required_argument, 0, 's'},
       {"ploidy", required_argument, 0, 'p'},
       {"distmat", required_argument, 0, 'd'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string genofile;
    string samplefile;
    string ploidyfile;
    bool write_distmat = false;
    string distmatfile;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "g:s:p:d:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'g':
                genofile = optarg;
                break;
            case 's':
                samplefile = optarg;
                break;
            case 'p':
                ploidyfile = optarg;
                break;
            case 'd':   
                write_distmat = true;
                distmatfile = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (genofile.length() == 0 || samplefile.length() == 0 || ploidyfile.length() == 0){
        fprintf(stderr, "ERROR: one or more required files was not provided\n");
        exit(1);
    }
    
    vector<string> samplenames;
    vector<int> ploidy;
    parse_samples(samplefile, samplenames);
    parse_ploidies(ploidyfile, ploidy);
    
    if (samplenames.size() != ploidy.size()){
        fprintf(stderr, "ERROR: different numbers of samples than ploidies: %ld %ld\n",
            samplenames.size(), ploidy.size());
        exit(1);
    }
    
    // Create distance matrix
    // Declare a distance matrix to use for building the UPGMA tree
    // Initialize to zero        
    vector<vector<float> > dist_mat;
    init_distmat(dist_mat, samplenames.size());
    
    vector<vector<float> > dist_mat_counts;
    init_distmat(dist_mat_counts, samplenames.size());
    
    // Read the input geno file
    
    fstream fin;
    fin.open(genofile.c_str(), fstream::in);
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", genofile.c_str());
        exit(1);
    }
    string line;
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        strip(line);
        
        if (line.length() > 0){
            if (line.length()  != samplenames.size()){
                fprintf(stderr, "Invalid number of samples in input file: %ld vs %ld\n",
                    line.length(), samplenames.size());
                exit(1);
            }
        
            for (int i = 0; i < line.length()-1; ++i){
                int alt_i = char2gt(line[i]);
                int ref_i = ploidy[i] - alt_i;
                
                if (alt_i == 9){
                    // Missing data
                    continue;
                }
                
                for (int j = i+1; j < line.length(); ++j){
                    // Difference is the total number of non-matching haplotypes.
                    int alt_j = char2gt(line[j]);
                    int ref_j = ploidy[j] - alt_j;
                    
                    if (alt_j == 9){
                        // Missing data
                        continue;
                    }
                    
                    int diff = alt_i - alt_j;
                    if (diff < 0){
                        diff = -diff;
                    }
                    dist_mat[i][j] += diff;
                    dist_mat_counts[i][j]++;
                    /*
                    int diff_ref = ref_i - ref_j;
                    if (diff_ref < 0){
                        diff_ref = -diff_ref;
                    }
                    int diff_alt = alt_i - alt_j;
                    if (diff_alt < 0){
                        diff_alt = -diff_alt;
                    }
                    dist_mat[i][j] += (diff_ref + diff_alt);
                    */
                }
    
            }
        }
    }
    
    // Normalize dist mat
    for (int i = 0; i < samplenames.size()-1; ++i){
        for (int j = i + 1; j < samplenames.size(); ++j){
            if (dist_mat_counts[i][j] > 0){
                dist_mat[i][j] /= dist_mat_counts[i][j];
            }
        }
    }
    
    // Now build the tree
    upgma(samplenames.size(), dist_mat, write_distmat, distmatfile, samplenames);
    

    return 0;
}

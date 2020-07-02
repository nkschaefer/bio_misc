// Defines treeNode class to represent pieces of trees for ARG inference.

#include <string>
#include <set>
#include <vector>
#include <stdio.h>
#include "treeNode.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <cmath>
#include <bitset>
#include <zlib.h>
#include <unordered_set>

using namespace std;

/**
 * Converts a bitset to a set of indices with value == 1
 */
set<unsigned int> bitset2set(const cladeset& b, const unsigned int num_haplotypes){
    set<unsigned int> s;
    // This is important: bitsets are stored in reverse order from the string
    // that created them. This means that haplotype indices are actually 
    // num_haplotypes - 1 - bit index.
    for (int i=0; i < num_haplotypes; i++){
        if (b[i]){
            s.insert(num_haplotypes-1-i);
        }
    }
    return s;
}

/**
 * Converts a set of indices to a bitset with those indices set to 1
 */
cladeset set2bitset(set<unsigned int>& b, const unsigned int num_haplotypes){
    cladeset s;
    for (set<unsigned int>::iterator it = b.begin(); it != b.end(); ++it){
        s.set(num_haplotypes-(*it)-1);
    }
    return s;
}

/**
 * Performs set intersection (bitset)
 */
cladeset set_int_bitset(const cladeset& set1, const cladeset& set2){
    return set1 & set2;
}

/**
 * Performs set difference (bitsets)
 */
cladeset set_diff_bitset(const cladeset& set1, const cladeset& set2){
    return set1 & ~set2;
}

/**
 * Computes the union of two sets (bitsets)
 */
cladeset set_union_bitset(const cladeset& set1, const cladeset& set2){
    return set1 | set2;
}
/**
 * Initialize a treeNode object.
 */
 
/**
 * Empty constructor
 */
treeNode::treeNode(){
    persistence = -1;
    support = 0;
    prop = true;
    visited = false;
    subtree_cache = false;
    
    recomb_left = false;
    recomb_right = false;
    leaving_clade_left = false;
    leaving_clade_right = false;
    
    top_level_prop = false;
    present_when_loaded = false;
    
    recomb_left_subtree = false;
    recomb_right_subtree = false;
    
    parent = NULL;
    cladeset leaves;
    cladeset leaves_below;
    num_haps = 0;
    dist = 1;
    
    keep_children = false;
    
    children.reserve(20);
}

/**
 * Initialize with leaves
 */
treeNode::treeNode(const unsigned int num_haplotypes, const cladeset& add_leaves){
    
    persistence = -1;
    support = 0;
    prop = true;
    visited = false;
    subtree_cache = false;
    recomb_left = false;
    recomb_right = false;
    
    leaving_clade_left = false;
    leaving_clade_right = false;
    
    top_level_prop = false;
    present_when_loaded = false;
    
    recomb_left_subtree = false;
    recomb_right_subtree = false;
    
    keep_children = false;
    
    parent = NULL;
 
    this->leaves = add_leaves;
    
    dist = 1;
    num_haps = num_haplotypes;
}

/**
 * Build bitsets once size is known
 */
void treeNode::set_haps(const int num_haplotypes){
    this->num_haps = num_haplotypes;
}

void treeNode::set_dist_above(float dist_down){
    dist_down += this->dist;
    
    this->dist_above = dist_down;
    
    float denom = this->dist_above + this->dist_below;
    if (denom == 0){
        denom = 1;
    }
    this->dist_norm = this->dist / denom;
    
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        (*child)->set_dist_above(dist_down);
    }
}

void treeNode::add_dist(float dist){
    
    float denom = (float)this->children.size();
    if (denom == 0){
        denom = 1.0;
    }
    dist = dist * (1.0/denom);
    
    this->dist += dist;
    this->dist_above += dist;
    this->dist_norm = this->dist / (this->dist_above + this->dist_below);
    
    // This should be a leaf node (or as close to one as exists), so don't propagate down.
    
    // Propagate this change up to the top.
    if (this->parent != NULL){
        treeNode* parent = this->parent;
        bool stop = false;
        while (!stop){
            dist = dist * (1.0 / (float)parent->children.size());
            parent->dist_below += dist;
            float denom2 = parent->dist_above + parent->dist_below;
            if (denom2 == 0){
                denom2 = 1.0;
            }
            parent->dist_norm = parent->dist / denom;
            if (parent->parent == NULL){
                stop = true;
            }
            else{
                parent = parent->parent;
            }
        }
    }
}

/**
 * Make current node the root of the tree
 */
void treeNode::setRoot(){
    
    this->parent = NULL;
    
    this->set_dist_above(this->dist);
}

/**
 * delete a treeNode object.
 */
treeNode::~treeNode(){
    // Remove all leaves
    /**
    this->leaves.clear();
    
    this->leaves_below.clear();
    
    this->children.clear();
    **/
    
    //if (!this->keep_children){
        this->leaves.reset();
        this->leaves_below.reset();
    
        for (vector<cladeset >::iterator pc = this->partner_clades_left.begin();
            pc != this->partner_clades_left.end();){
            pc->reset();
            this->partner_clades_left.erase(pc);
        }
        for (vector<cladeset >::iterator pc = this->partner_clades_right.begin();
            pc != this->partner_clades_right.end();){
            pc->reset();
            this->partner_clades_right.erase(pc);
        }
    
        this->delete_children();
        this->children.clear();
    //}
}

/**
 * Deletes all treeNodes in the tree below a given node, including this node.
 */
void treeNode::free_recursive(){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ){
        (*child)->free_recursive();
        this->children.erase(child);
    }
    delete (this);
}

/**
 * Deletes all treeNodes in the tree below a given node, not including this
 * node.
 */
void treeNode::delete_children(){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end();){
        (*child)->free_recursive();
        //(*child)->delete_children();
        this->children.erase(child);
        //delete *child;
    }
}

/**
 * Copy constructor.
 */
treeNode::treeNode(const treeNode &t){
    
    this->support = t.support;
    this->visited = t.visited;    
    this->num_haps = t.num_haps;

    // Force regeneration of cached leaves in subtree.
    this->subtree_cache = false;
    //set<string> leaves_below;
    
    if (t.num_haps > 0){
        
        this->leaves = t.leaves;
    }
    
    else{
       
    }
    this->persistence = t.persistence;
    this->keep_children = false;
    this->top_level_prop = t.top_level_prop;
    this->prop = t.prop;
    this->dist = t.dist;
    this->parent = t.parent;
    this->recomb_left = t.recomb_left;
    this->recomb_right = t.recomb_right;
    
    this->leaving_clade_left = t.leaving_clade_left;
    this->leaving_clade_right = t.leaving_clade_right;
    
    this->present_when_loaded = t.present_when_loaded;
    
    this->recomb_left_subtree = t.recomb_left_subtree;
    this->recomb_right_subtree = t.recomb_right_subtree;
    
    this->partner_clades_left = t.partner_clades_left;
    this->partner_clades_right = t.partner_clades_right;
    
    this->name = t.name;
    
    // Copy all children as well.
    vector<treeNode*> children;
    for (int i=0; i < t.children.size(); ++i){
        this->children.push_back(new treeNode((*t.children[i])));
    }
}

/**
 * operator= function
 * Given a treeNode and another treeNode, modifies this treeNode so it exactly
 * resembles the given treeNode
 */
treeNode& treeNode::operator=(const treeNode &t){
    // Check for self assignment
    if (this == &t){
        return *this;
    }
    else{
        this->persistence = t.persistence;
        this->keep_children = false;
        this->support = t.support;
        this->visited = t.visited;
        this->dist = t.dist;
        this->prop = t.prop;
        this->top_level_prop = t.top_level_prop;
        
        this->recomb_left_subtree = t.recomb_left_subtree;
        this->recomb_right_subtree = t.recomb_right_subtree;
        
        if (t.num_haps == 0){
            // If this happens, we'll need to set num_haps at a later step. Otherwise,
            // the bitsets will throw errors.
            
            this->num_haps = 0;
        }
        else if (this->num_haps == t.num_haps){
            // No need to re-declare; just re-set.
            this->leaves.reset();
            this->leaves = t.leaves;
            this->leaves_below.reset();
        }
        else{
            
            this->num_haps = t.num_haps;
        }
        
        // Force regeneration of cached leaves in subtree.
        this->subtree_cache = false;
        
        this->recomb_left = t.recomb_left;
        this->recomb_right = t.recomb_right;
        this->leaving_clade_left = t.leaving_clade_left;
        this->leaving_clade_right = t.leaving_clade_right;
        
        this->present_when_loaded = t.present_when_loaded;
        
        this->partner_clades_left = t.partner_clades_left;
        this->partner_clades_right = t.partner_clades_right;
        
        this->name = t.name;
        
        // Erase children
        this->children.clear();
        
        // Copy all children over
        vector<treeNode*> children;
        for (int i=0; i < t.children.size(); ++i){
            this->children.push_back(new treeNode((*t.children[i])));
        }
        
        return *this;
    }
}

/**
 * Given a node representing the root of a tree, recursively visits all that node's 
 * children and sets parent pointers.
 */
void add_parents(treeNode& t){
    for (vector<treeNode*>::iterator child = t.children.begin(); child != t.children.end();
        ++child){
        (*child)->parent = &t;
        add_parents(*(*child));
    }
}

/**
 * Removes cached information about all leaves in this node's subtree.
 */
void treeNode::clear_cache(){
    //this->leaves_below = set<string>();
    this->leaves_below.reset();
    this->subtree_cache = false;
}

/**
 * Adds a child to the tree below this node; notifies that child that this
 * node is now its parent.
 *
 * http://codereview.stackexchange.com/questions/47395/tree-template-class-implementation-for-c
 */
void treeNode::addChild(cladeset& leaves, float brlen){
    treeNode* child = new treeNode(this->num_haps, leaves);
    child->dist = brlen;
    children.push_back(child);
    child->parent = this;
    
    this->leaves = set_diff_bitset(this->leaves, leaves);
    this->clear_cache();
}

void treeNode::addLeaf(const unsigned int leaf){
    this->leaves.set(this->num_haps-1-leaf);
    this->clear_cache();
}

/**
 * Return a Newick-format string representation of this treeNode.
 * 
 * Arguments:
 *  recomb -- true or false, whether or not this node's status as marked by
 *  recombination (or not) in either direction should be marked in the Newick string.
 *  use_hapnames -- true or false -- should haplotype indices be converted into names
 *      using the provided (possibly empty) haplotype index -> name mapping?
 *  hapnames -- a mapping of haplotype index -> string name to use in the tree
 */
string treeNode::newick(bool recomb, bool use_hapnames, vector<std::string> hapnames){
    string recomb_left_str = "";
    
    if (recomb && this->leaving_clade_left){
        recomb_left_str = "{";
    }
    
    if (recomb && this->recomb_left){
        recomb_left_str += "[";
    }
    string recomb_right_str = "";
    
    if (recomb && this->leaving_clade_right){
        recomb_right_str = "}";
    }
    
    if (recomb && this->recomb_right){
        recomb_right_str += "]";
    }
    string wrapperL = "";
    string wrapperR = "";
    string root = "";
    if (this->parent == NULL){
        wrapperL = "(";
        wrapperR = ")";
        
        char rootbuf[15];
        sprintf(rootbuf, ":%5.5f;", this->dist_norm);
        root += string(rootbuf);
        
        //root = ";";
    }

    if (this->children.size() == 0){
        // Base case: this is just leaves.
        string newickstr = "";
        //if (false){
        if (leaves.count() == 1){
            set<unsigned int> leaves_set = bitset2set(leaves, this->num_haps);
            if (use_hapnames){
                newickstr = hapnames[*leaves_set.begin()];
            }
            else{
                char leaf_str[20];
                sprintf(leaf_str, "%d", *leaves_set.begin());
                newickstr = string(leaf_str);
            }
        }
        else if (leaves.count() > 0){
            set<unsigned int> leaves_set = bitset2set(leaves, this->num_haps);
            set<unsigned int>::iterator lastit = leaves_set.end();
            lastit--;
            for (set<unsigned int>::iterator leaf = leaves_set.begin(); leaf != leaves_set.end(); ++leaf){
                if (use_hapnames){
                    string leaf_str = hapnames[*leaf] + ":0";
                    newickstr += leaf_str;
                }
                else{
                    char leaf_str[20];
                    sprintf(leaf_str, "%d:0", (*leaf));
                    newickstr += string(leaf_str);
                }
                //newickstr += (*leaf) + ":0";
                if (leaf != lastit){
                    newickstr += ",";
                }
            }
        }
        return (wrapperL + recomb_left_str + newickstr + recomb_right_str + wrapperR + root);
    }
    else{
        
        string children_rendered = "";
        vector<treeNode*>::iterator lastit = this->children.end();
        lastit--;
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end(); ++child){
            // Note: adding distances here (only when the method is called on a 
            // child node) means that no distance will be given to the root node.
            // That's okay, since the concept of distance doesn't make sense
            // for it.
            string child_newick = (*child)->newick(recomb, use_hapnames, hapnames);
            
            if (this->name.length() == 0){
                char buf[child_newick.size() + 10 + 3];
                
                //if (false){
                if ((*child)->leaves.count() == 1 && (*child)->children.size() == 0){
                    // Child is a leaf node.
                    sprintf(buf, "%s:%5.5f", child_newick.c_str(), abs((*child)->dist_norm));
                }
                else{
                    sprintf(buf, "(%s):%5.5f", child_newick.c_str(), abs((*child)->dist_norm));
                }
                children_rendered += (string) buf;
            }
            else{
                char buf[child_newick.size() + this->name.length() + 10 + 3 + 1];
                
                sprintf(buf, "(%s)%s:%5.5f", child_newick.c_str(), this->name.c_str(), abs((*child)->dist_norm));
            
                children_rendered += (string) buf;
            }
            if (child != lastit || this->leaves.count() > 0){
                children_rendered += ",";
            }
        }
        
        string leaves_rendered = "";
        if (this->leaves.count() > 0){
            set<unsigned int> leaves_set = bitset2set(this->leaves, this->num_haps);
            set<unsigned int>::iterator lastleaf = leaves_set.end();
            //set<string>::iterator lastleaf = this->leaves.end();
            lastleaf--;
            for (set<unsigned int>::iterator leaf = leaves_set.begin(); leaf != leaves_set.end(); ++leaf){
                if (use_hapnames){
                    string leaf_str = hapnames[*leaf];
                    leaves_rendered += leaf_str + ":0";
                }
                else{
                    char leaf_str[12];
                    sprintf(leaf_str, "%d:0", (*leaf));
                    leaves_rendered += string(leaf_str);
                }
                //leaves_rendered += (*leaf) + ":0";
                if (leaf != lastleaf){
                    leaves_rendered += ",";
                }
            }
        }
          
        return (wrapperL + recomb_left_str + children_rendered + leaves_rendered + 
            recomb_right_str + wrapperR + root);
    }
}

/**
 * Return a vector of pointers to all treeNodes representing lowest-level nodes
 * (no children) below this node in the tree.
 */
set<treeNode*> treeNode::get_lowest(){
    set<treeNode*> lowest_nodes;
    if (children.size() == 0){
        lowest_nodes.insert(this);
    }
    else{
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end(); ++child){
            set<treeNode*> child_lowest_nodes = (*child)->get_lowest();
            for (set<treeNode*>::iterator child_child = child_lowest_nodes.begin();
                child_child != child_lowest_nodes.end(); ++child_child){
                lowest_nodes.insert(*child_child);
            }
        }
    }
    return lowest_nodes;
}

void treeNode::clear_cache_above(){
    this->clear_cache();
    if (this->parent == NULL){
        return;
    }
    else{
        this->parent->clear_cache_above();
    }
}

/**
 * Return a set of all leaf names (strings) in this node's entire subtree.
 */
cladeset treeNode::subtree_leaves(){
    //if (this->leaves_below.size() == 0){
    //    this->leaves_below = cladeset(this->num_haps);
    //}
    if (subtree_cache){
        // We have a stored version of this set; no need to re-generate.
        return this->leaves_below;
    }
    else{
        // Need to re-generate.
        this->leaves_below = this->leaves;
        for (vector<treeNode*>::iterator child = this->children.begin();
            child != this->children.end(); ++child){
            cladeset child_subtree = (*child)->subtree_leaves();
            this->leaves_below = set_union_bitset(this->leaves_below, child_subtree);
        }
        this->subtree_cache = true;
        return this->leaves_below;
    }
}


/**
 * Called by parse_newick().
 * Finds all Newick strings representing subtrees at the current depth of the
 * Newick string and returns them as a vector of strings. Each can then
 * be parsed independently by parse_newick().
 */
void parse_subtree(vector<string>& subtrees, string& newick){
    //vector<string> subtrees;
    if (newick.size() == 0){
        return;
        //return subtrees;
    }
    
    // Store indices of subtrees.
    vector<pair<int, int> > subtree_inds;
    
    int parenCount = 0;
    int last_subtree_start = 0;
    int last_subtree_end = 0;
    
    for (int chr_index=0; chr_index < newick.size(); ++chr_index){
        if (newick[chr_index] == '('){
            parenCount++;
        }
        else if (newick[chr_index] == ')'){
            parenCount--;
            
            if (parenCount == 0){
                // Two things can happen here: either the string ends,
                // or we find a comma to indicate that this is the
                // start of another subtree at the same level.
                int subtree_end = chr_index+1;
                
                for (int chr_index2 = chr_index+1; chr_index2 < newick.size(); ++chr_index2){
                    if (newick[chr_index2] == ','){
                        subtree_end = chr_index2;
                        break;
                    }
                    else if (chr_index2 == newick.size()-1){
                        // We've hit the end of the tree without finding
                        // a comma.
                        subtree_end = newick.size();
                    }
                }
                
                pair<int, int> inds = make_pair(last_subtree_start, subtree_end);
                subtree_inds.push_back(inds);
                
                // Start next subtree at this index so we can avoid hitting
                // the comma that separates subtrees
                last_subtree_start = subtree_end + 1;
            }
        }
    }
    // Note: string must end with closed parenthesis, so the final subtree
    // indices should already be added to the list.
    for (vector<pair<int, int> >::iterator subtree_ind = subtree_inds.begin();
        subtree_ind != subtree_inds.end(); ++subtree_ind){
        char subtree[subtree_ind->second - subtree_ind->first + 1];
        memcpy(subtree, newick.c_str() + subtree_ind->first, subtree_ind->second - subtree_ind->first);
        subtree[subtree_ind->second - subtree_ind->first] = '\0';
        string subtree_str = string(subtree);
        subtrees.push_back(subtree_str);   
    }
        
    //return subtrees;
}



string extract_regex_match(string subject, const std::sub_match<const char*>& fullmatch, 
    const std::sub_match<const char*>& groupmatch){
    
    long int startpos = groupmatch.first - fullmatch.first;
    long int cpylen = groupmatch.second - groupmatch.first;
    
    char extracted[cpylen+1];
    memcpy(extracted, subject.c_str() + startpos, cpylen);
    extracted[cpylen] = '\0';
    return (string) extracted;
    
}


/**
 * Called by parse_newick().
 * Given a string representing a leaf, parses the name and branch length of the leaf.
 * Determines whether or not the leaf deserves its own node and either creates a node
 * and adds it as a child of the given node, or adds the leaf directly to the node's
 * set of leaves (modifies parent by reference)
 */
bool parse_leaf(string& leaf, string& leaf_add, float& brlen){
    // Regular expression for parsing leaves (and associated branch lengths)
    static regex leafmatch(R"((.+):([0-9\.\-]+))", regex::extended);
    
    // Return true = add leaf as child
    // or false = add leaf as leaf
    cmatch match;
    if (regex_match(leaf.c_str(), match, leafmatch)){
        string leafname = extract_regex_match(leaf, match[0], match[1]);
        string leafdist = extract_regex_match(leaf, match[0], match[2]);
        
        float dist = abs(atof(leafdist.c_str()));
       
        if (dist > 0){
            if (leafname.length() > 0){
                // Create a new node.
                leaf_add = leafname;
                brlen = dist;
                return true;
            }
        }
        else{
            leaf_add = leafname;
            return false;
        }
    }
    else{
        // No branch length provided.
        leaf_add = leaf;
        return false;
    }
    return false;
}

/**
 * (Recursive method)
 * Given a string in Newick format, parses it as a tree. Note that
 * trees nested within this string are also trees in Newick format.
 *
 * leaf_branches = allow branch lengths for leaves (default behavior is to
 *    create a separate node for each leaf that has a nonzero branch
 *    length. If False, all leaves will be treated as branch length zero
 *    and will not be given their own nodes.
**/
void parse_newick(treeNode& root, string newick, const unsigned int num_haplotypes){
    
    root.set_haps(num_haplotypes);
    root.dist = 0.0;
    
    bool is_root = false;
    if (newick[newick.size()-1] == ';'){
        is_root = true;
        newick.erase(newick.end()-1);
    }

    // Precompile all regular expressions
    // NOTE: we use Boost library regexes here, since they're Perl-style.
    try{
        // Regular expression to match whole Newick-format tree nodes (with final 
        // semicolon removed)
        //static regex nodematch(R"(\(([a-zA-Z0-9:\.,\(\)\[\]\{\}\-_]+)\)([a-zA-Z0-9\-_]+)?(:[0-9\.]+)?)", regex::extended);
        static regex nodematch(R"(\(([a-zA-Z0-9:.,()\-_]+)\)([a-zA-Z0-9\-_]+)?(:[0-9\.e\-]+)?)", regex::extended);

        // Get all root-level information out of the tree.
        // Note: we have to look for brackets inside strings, on the left
        // or right side. These indicate clades being blocked off by
        // recombination.
    
        cmatch match;
        if (regex_match(newick.c_str(), match, nodematch)){
    
            // Innertree will be whatever occurs inside the outermost set of 
            // parenthesis.
            string innertree = extract_regex_match(newick.c_str(), match[0], match[1]);
        
            // Look for recombination marks.
        
            // Look for node name
            if (match[2].second != match[2].first){
                //root.name = extract_regex_match(newick.c_str(), match[0], match[2]);
            }
        
            // Look for branch length
            if (match[3].second != match[3].first){
                string dist_str = extract_regex_match(newick.c_str(), match[0], match[3]);
                if (dist_str[0] == ':'){
                    dist_str = dist_str.substr(1, dist_str.length()-1);
                }
                root.dist = abs(atof(dist_str.c_str()));
            }
        
            // Find all leaves attached to this node.
        
            // First, check for leaves on the left.
        
            const char *left_paren_ptr = strchr(innertree.c_str(), '(');
            if (left_paren_ptr && left_paren_ptr != innertree.c_str()){
                int left_paren_index = left_paren_ptr - innertree.c_str();
                char leafStr[left_paren_index+1];
                strncpy(leafStr, innertree.c_str(), left_paren_index);
                leafStr[left_paren_index] = '\0';
                if (leafStr[0] == ','){
                    int trunc_index = strlen(leafStr)-1;
                    memmove(leafStr, leafStr+1, strlen(leafStr)-1);
                    leafStr[trunc_index] = '\0';
                }
                if (leafStr[strlen(leafStr)-1] == ','){
                    leafStr[strlen(leafStr)-1] = '\0';
                }
            
                string leafStr_str(leafStr);
                istringstream tokenizer(leafStr_str);
                string leaf;
            
                vector<string> leaves_str;
                while(std::getline(tokenizer, leaf, ',')){
                    leaves_str.push_back(leaf);
                }
                
                for (vector<string>::iterator ls = leaves_str.begin(); ls != leaves_str.end(); ++ls){
                    string leafstr = "";
                    float brlen = 0;
                    bool add_leaf = parse_leaf(*ls, leafstr, brlen);
                    unsigned int leaf_index = atoi(leafstr.c_str());
                    if (add_leaf){
                        if (false){
                        //if (leaves_str.size() == 1){
                            root.leaves.set(root.num_haps-1-leaf_index);
                            root.dist = brlen;
                        }
                        else{
                            cladeset leaves;
                            leaves.set(root.num_haps-1-leaf_index);
                            root.addChild(leaves, brlen);
                        }
                    }
                    else{
                        root.leaves.set(root.num_haps-1-leaf_index);
                    }
                }
                
                // Pare tree down to part after leaves
                innertree.erase(0, left_paren_index);
            
            }
        
            // Now check for leaves on the right.
            const char *right_paren_ptr = strrchr(innertree.c_str(), ')');
            if (right_paren_ptr && (right_paren_ptr - innertree.c_str()) != strlen(innertree.c_str())){
                int right_paren_index = right_paren_ptr - innertree.c_str();
                char leafStr[strlen(innertree.c_str())-right_paren_index+1];
                memcpy(leafStr, innertree.c_str()+right_paren_index+1, strlen(innertree.c_str())-right_paren_index);
                leafStr[strlen(innertree.c_str())-right_paren_index] = '\0';
            
                // Figure out whether the first "leaf" is actually a branch length.
                if (leafStr[0] == ':'){
                    int brlen_end_index;
                    for (int i = 0; i < strlen(leafStr); ++i){
                        brlen_end_index = i;
                        if (leafStr[i] == ','){
                            break;
                        }
                    }
                    // Include branch length in inner tree.
                    innertree.erase(right_paren_index + brlen_end_index, innertree.length() - right_paren_index - brlen_end_index - 1);
                
                    if (brlen_end_index < strlen(leafStr)-1){
                        // Chop off the branch length portion from leafStr.
                        int trunc_index = strlen(leafStr) - brlen_end_index - 1;
                        memmove(leafStr, leafStr + brlen_end_index + 1, strlen(leafStr) - brlen_end_index-1);
                        leafStr[trunc_index] = '\0';
                    }
                    else{
                        // leafStr is just a branch length 
                        // do nothing -- no leaves to process.
                        leafStr[0] = '\0';
                    }
                }
                else{
                    innertree.erase(right_paren_index+1, innertree.length()-right_paren_index+1);
                }
                // Only proceed if the leaf string is not empty
                if (strlen(leafStr) > 0){
                    string leafStr_str(leafStr);
            
                    istringstream tokenizer(leafStr_str);
                    string leaf;
                
                    while(std::getline(tokenizer, leaf, ',')){
                        string leafstr = "";
                        float brlen = 0;
                        bool add_leaf = parse_leaf(leaf, leafstr, brlen);
                    
                        unsigned int leaf_index = atoi(leafstr.c_str());
                        if (add_leaf){
                            cladeset leaves;
                            leaves.set(root.num_haps-1-leaf_index);
                            root.addChild(leaves, brlen);
                        }
                        else{
                            root.leaves.set(root.num_haps-1-leaf_index);
                        }
                    }
                }
            
                // We now have a leafless tree. Look for all clades contained within
                // and parse each.
            
                // We now have a smaller Newick-format tree or set of trees.
                // parse_subtree will look to see if it's one or multiple trees
                // at this level and parse each.
                vector<string> subtrees;
                parse_subtree(subtrees, innertree);
                
                for (vector<string>::iterator subtree_it = subtrees.begin();
                    subtree_it != subtrees.end(); ++subtree_it){
                    treeNode* child = new treeNode;
                    child->parent = &root;
                    child->set_haps(root.num_haps);
                    parse_newick(*child, (*subtree_it), num_haplotypes);
                    // NOTE: don't touch this. You have to make a copy or else
                    // you'll store a reference to something that gets deleted later on.
                    //root.children.push_back(new treeNode(child));
                
                    root.children.push_back(child);
                }
            
            }
            else{
                // Base case: this tree was just leaves.
                istringstream tokenizer(innertree);
                string leaf;
                while(std::getline(tokenizer, leaf, ',')){
                    string leafstr = "";
                    float brlen = 0;
                    bool add_leaf = parse_leaf(leaf, leafstr, brlen);
                
                    unsigned int leaf_index = atoi(leafstr.c_str());
                    if (add_leaf){
                        cladeset leaves;
                        leaves.set(root.num_haps-1-leaf_index);
                        root.addChild(leaves, brlen);
                    }
                    else{
                        root.leaves.set(root.num_haps-1-leaf_index);
                    }
                }
            }
        }
        else{
            // Not valid Newick format.
            fprintf(stderr, "ERROR parsing tree:\n");
            fprintf(stderr, "%s\n", newick.c_str());
            exit(1);
        }
    } catch (regex_error& e){
        if (e.code() == regex_constants::error_collate){
            fprintf(stderr, "regex error_collate\n");
        }
        else if (e.code() == regex_constants::error_ctype){
            fprintf(stderr, "regex error_ctype\n");
        }
        else if (e.code() == regex_constants::error_escape){
            fprintf(stderr, "regex error_escape\n");
        }
        else if (e.code() == regex_constants::error_backref){
            fprintf(stderr, "regex error_backref\n");
        }
        else if (e.code() == regex_constants::error_brack){
            fprintf(stderr, "regex error_brack\n");
        }
        else if (e.code() == regex_constants::error_paren){
            fprintf(stderr, "regex error_paren\n");
        }
        else if (e.code() == regex_constants::error_brace){
            fprintf(stderr, "regex error_brace\n");
        }
        else if (e.code() == regex_constants::error_badbrace){
            fprintf(stderr, "regex error_badbrace\n");
        }
        else if (e.code() == regex_constants::error_range){
            fprintf(stderr, "regex error_range\n");
        }
        else if (e.code() == regex_constants::error_space){
            fprintf(stderr, "regex error_space\n");
        }
        else if (e.code() == regex_constants::error_badrepeat){
            fprintf(stderr, "regex error_badrepeat\n");
        }
        else if (e.code() == regex_constants::error_complexity){
            fprintf(stderr, "regex error_complexity\n");
        }
        else if (e.code() == regex_constants::error_stack){
            fprintf(stderr, "regex error_stack\n");
        }
        exit(1);
    }
    

    //return root;    
}


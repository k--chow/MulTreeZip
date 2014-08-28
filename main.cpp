/*
 * A project to compress/decompress weighted or unweight  multi-labeled trees (mul-trees)
 * the mul-trees can have overlapping leaf multisets
 * 
 * File:   main.cpp
 * Author: ruchi Chaudhary
 *
 * Created on May 7, 2013, 4:04 PM
 */

#include <stdlib.h>
#include <bitset>
#include "common.h"
#include "argument.h"
#include "gauge.h"
#include "global.h"
#include "parsing.h"
#include "tree.h"
#include "tree_IO.h"
#include "tree_traversal.h"
#include "tree_subtree_info.h"
#include "tree_name_map.h"
#include "hash.h"
#include "compressfunc.h"
#include "decompressfunc.h"
#include "isomor.c"

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>

#ifndef NOHASH
#include <boost/unordered_set.hpp>
#endif


typedef boost::unordered_map<unsigned int,int> gid2ctype;

/*
 * 
 */
int main(int ac, char* av[]) {
    {
        std::ostringstream os; os << "command:";
        for (int i = 0; i < ac; i++) os << ' ' << av[i];
        MSG(os.str());
    }
    std::string in_file;    
    std::string out_file;    
    bool compress = false;
    bool decompress = false;
    bool uniq = false;
    unsigned int seed = std::time(0);
    
    {
        Argument a; a.add(ac, av);
        // help
        if (a.existArg2("-h","--help")) {
            MSG("options:");
            MSG("  -i [ --input ] arg      input trees (.newick for --compress and .mtz for decompress)");
            MSG("  -o [ --output ] arg     output file (.mtz for --compress and .newick for decompress)");
            MSG("       -c         compress the input .newick file into .mtz file");
            MSG("       -d         decompress the input .mtz file into .newick file");
            MSG("       -u         only unique");
            MSG("       --seed arg         random generator seed");
            MSG("  -h [ --help ]           produce help message");
            MSG("");
            MSG("example:");
            MSG("  " << av[0] << " -i infile.newick --compress -o outfile.mtz");
            exit(0);
        }

        // input trees
        if (a.existArgVal2("-i", "--input", in_file)) MSG("input file: " << in_file) else MSG("using standard input");
        
        // output file
        if (a.existArgVal2("-o", "--output", out_file)) MSG("output file: " << out_file);
        
        // compress the input newick file
        compress = a.existArg("-c");

        // compress uniques and list all 
        uniq = a.existArg("-u");
        
        // decompress the input .mtz file
        decompress = a.existArg("-d");

        if((compress && decompress) || (!compress && !decompress))
            ERROR_exit("Wrong argument: choose either compress or decompress");
        
        // random seed
        a.existArgVal("--seed", seed);
        MSG("seed: " << seed);
        //aw::rng.seed(static_cast<unsigned int>(seed));
        // unknown arguments?
        a.unusedArgsError();
    }
    //------------------------------------------------------------------------------------------------------------------


    int t3 = clock();    

    // ------------------If for Compression Begins-------------------------------------------------------------------------
    if(compress) {
        // read input trees
        std::vector<aw::Tree> g_trees;
        std::vector<aw::idx2name> g_taxa;
        std::vector<std::vector<long double> > wght;        
        HashRFMap vvec_hashrf;
        
        {
            const std::string filename = in_file;
            { // read trees
                std::ifstream ifs;
                std::istream &is = filename.empty() ? std::cin : ifs;
                if (!filename.empty()) {
                    ifs.open(filename.c_str());
                    if (!ifs) ERROR_exit("cannot read file '" << filename << "'");
                }

                MSG_nonewline("Reading input trees: ");
                aw::gauge_exp g; aw::gauge_init(&g);
                for (;;) {
                    aw::Tree t;
                    aw::idx2name t_names;
                    std::vector<long double> t_weight;
                    if (!aw::stream2tree(is, t, t_names,t_weight)) break;
                    g_taxa.push_back(t_names);
                    g_trees.push_back(t);
                    if(WEIGHTED) {
                        if(t_weight.empty()) WEIGHTED = false;
                        else  wght.push_back(t_weight);
                    }
                    aw::gauge_inc(&g);
                }
                aw::gauge_end(&g);
            }
            MSG("Input trees: " << g_trees.size());
            if (g_trees.empty()) ERROR_exit("No input trees found in file '" << filename << "'");
        }

        if(WEIGHTED)
            MSG("Weighted trees...")
        else
            MSG("Unweighted trees...");

        // map taxa labels
        aw::TaxaMap taxamap;  //to store all taxon and global id
        std::vector<aw::TreetaxaMap> g_nmaps; //for mapping taxamap and idx2name (of a tree)
        {
            g_nmaps.resize(g_taxa.size());
            for (unsigned int i=0,iEE=g_taxa.size(); i<iEE; ++i) {
                aw::idx2name &n = g_taxa[i];
                taxamap.insert(n);
                g_nmaps[i].create(n,taxamap);
            }
            MSG("Taxa: " << taxamap.size());
        }

//        std::ofstream ouput_fs;
//        if (!out_file.empty()) {
//            ouput_fs.open(out_file.c_str());
//            if (!ouput_fs) ERROR_exit("cannot write file '" << out_file << "'");
//        }
//        std::ostream &output = out_file.empty() ? std::cout : ouput_fs;
        
        //HOW MANY UNIQUE TREES --------------------------------------------------------------------------------------------------------------------------------------------
        unq_tr.push_back(0);  //g_trees[0]

        //for each g_trees if similar to a unique tree
        sim2unq.push_back(0);  //g_trees[0] is similar to uniq_tnum[0]

        MSG_nonewline("Processing trees...");

        //aw::tree2newick(output,g_trees[0],g_taxa[0]); output << std::endl;
        
        if(!WEIGHTED && uniq) {
            for(unsigned int i=1,iEE=g_trees.size(); i<iEE; ++i) {
                bool flag = true;
                //std::cout<<"UNS:"<<unq_num.size();
                for(unsigned int j=0,jEE=unq_tr.size(); j<jEE; ++j) {
                    unsigned int k = unq_tr[j];
                    //std::cout<<"ik"<<i<<" "<<k<<".";
                    int is = isomorphism(g_trees[i],g_trees[k],g_nmaps[i],g_nmaps[k]);
                    if(is == 0) {  //isomorphism
                        //std::cout<<"<SIM>";
                        sim2unq.push_back(j);
                        flag = false;
                        break;
                    }
                }
                if(flag) { //no similarity in the unique list
                    //aw::tree2newick(output,g_trees[i],g_taxa[i]); output << std::endl;
                    unq_tr.push_back(i);
                    sim2unq.push_back(unq_tr.size()-1);
                }                
            }
        }
        else{
            for(unsigned int i=1,iEE=g_trees.size(); i<iEE; ++i) {
                unq_tr.push_back(i);
                sim2unq.push_back(unq_tr.size()-1);
            }
        }


        MSG("\nUnique Trees: "<<unq_tr.size());

        //part for unique trees ------------------------------------
//        std::ofstream ouput_fs;
//        if (!out_file.empty()) {
//            ouput_fs.open(out_file.c_str());
//            if (!ouput_fs) ERROR_exit("cannot write file '" << out_file << "'");
//        }
//        std::ostream &output = out_file.empty() ? std::cout : ouput_fs;
//
//        output << "[Unique Trees = " <<unq_tr.size()<<"]"<<std::endl;
//        for(unsigned int i=0,iEE=unq_tr.size(); i<iEE; ++i) {
//            unsigned int k = unq_tr[i];
//            aw::tree2newick(output,g_trees[k],g_taxa[k]); output << std::endl;
//
//        }
//        return(0);
        //------------------------------------------------------------
        

        //WE WILL COMPRESS ONLY UNIQUE TREES--------------------------------------------------------------------------------------------------------------------------------

        //counting maximum number of copies of a gene node in a genetree        
        gid2ctype gid2cnt;    //<gid,count>
        {
            for(unsigned int i=0,iKK=taxamap.size(); i<iKK; ++i) {
                gid2cnt[i]=0;
                order.push_back(taxamap.taxon(i));
                //std::cout<<" "<<i<<" - "<<order[i];
                lab2order[order[i]] = i;
            }
            for (unsigned int ii=0,iiEE=unq_tr.size(); ii<iiEE; ++ii) {
                unsigned int i = unq_tr[ii];
                bool mulTree = false;
                TREE_FOREACHLEAF(v,g_trees[i]) {
                    unsigned int gid_v = g_nmaps[i].gid(v);
                    int gid_cnt = g_nmaps[i].ids_count(gid_v);
                    if(gid_cnt>1) mulTree = true;
                    if(gid2cnt[gid_v] < gid_cnt) gid2cnt[gid_v] = gid_cnt; }
                if(mulTree) all_singly = false;
                if(!mulTree) all_mul = false;
                multi_labeled.push_back(mulTree);
            }           
        }
      
        unsigned int bits = 0;        
        {
             BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
                std::string name = taxamap.taxon(w.first);
                int copies = w.second;

                std::string bin_string = decimal2binary(copies);
                int len_bin = bin_string.length();
                                             
                bits +=len_bin;                
                guide[name] = bits-1;   
            }            
        }

        //updating global variables
        NUM_TREES = g_trees.size();
        NUM_TAXA = taxamap.size(); // number of taxa
        BITSET = bits;       // number of taxa and their multiplicities

        std::cout<<"\nBITSET: "<<BITSET;

        initialize_hashtable(vvec_hashrf, M1, M2); //initialize contents of hashtable

        std::cout<<"\nM1 & M2:"<<M1<<" "<<M2;

        BIG_M1 = M1+1;
        
        //compress each tree
        {
            for (int ii=0,iiEE=unq_tr.size(); ii<iiEE; ++ii) {
                int i = unq_tr[ii];
                aw::Tree g_tree = g_trees[i];
                
                //for easy parent-child relationship in a g_tree
                aw::SubtreeParent<aw::Tree>  g_parent;                             
                g_parent.create(g_tree);

                boost::dynamic_bitset<> * bs_vec[g_tree.node_size()];
                dfs_compute_hash_w_bitstring(g_tree.root, g_tree, g_parent, g_taxa[i], vvec_hashrf, i, bs_vec);

                boost::unordered_map<unsigned int, unsigned long long> ord;
                for (aw::Tree::iterator_dfs m=g_tree.begin_dfs(),mEE=g_tree.end_dfs(); m!=mEE; ++m) {
                    if(m.idx == g_tree.root) continue;
                    unsigned long long par_ord, par_hv1, par_hv2;

                    switch (m.direction) {
                        case aw::PREORDER: {
                            if(m.parent==g_tree.root) {
                                par_ord = 0; par_hv1 = BIG_M1; par_hv2 = 0;
                            }
                            else {
                                par_ord = ord[m.parent]; par_hv1 = g_tree.return_hv1(m.parent); par_hv2 = g_tree.return_hv2(m.parent);
                            }
                            long double m_weight = 0;
                            if(WEIGHTED)  m_weight = wght[i][m.idx];
                            int s_idx = sim2unq[i];

                            //we will store those leaves that have their parents as root
                            bool par_root = false;
                            if(m.parent == g_tree.root)
                                par_root = true;
                            //std::cout<<"\n*"<<m.idx<<" "<<par_root;
                            unsigned long long clust_column =  vvec_hashrf.hashing_bs(s_idx,g_tree.is_leaf(m.idx),par_root,g_tree.return_hv1(m.idx),g_tree.return_hv2(m.idx),par_hv1,par_hv2,par_ord,m_weight,bs_vec[m.idx]);
                            ord[m.idx] = clust_column;
                        }
                    }                    
                }  
            }
            //store vvec_hashrf in a file
            compressf(out_file, vvec_hashrf);
        }       

        MSG("\nCompression successful!!!");
    }  // ---------- If of Compression Ends -----------------------------------------------------------------

    if(decompress) {
        std::ofstream ouput_fs;
        if (!out_file.empty()) {
            ouput_fs.open(out_file.c_str());
            if (!ouput_fs) ERROR_exit("cannot write file '" << out_file << "'");
        }
        std::ostream &output = out_file.empty() ? std::cout : ouput_fs;
        
        std::vector<aw::Tree> trees;
        std::vector<aw::idx2name> names;

        MSG("\nDecompressing...")
        decompression(in_file,trees,names);

        //writing trees in the file...
        for(int i=0, iEE=trees.size(); i<iEE; ++i) {
            aw::tree2newick(output,trees[i],names[i]); 
            output << std::endl;
        }
        
        MSG("\nDecompression successful!!!")
    } //-----------if for Decompression Ends here

    //for timing 
    int t4 = clock();
    {   
        long ttime = ((long)(t4-t3))/CLOCKS_PER_SEC;
        int d,h,m,s;
        util::convertTime(ttime,d,h,m,s);
        MSG_nonewline("\nTotal elapsed time: ");
        if(d!=0) MSG_nonewline(d<<"d ");
        if(h!=0) MSG_nonewline(h<<"h ");
        if(m!=0) MSG_nonewline(m<<"m ");
        MSG_nonewline(s<<"s ");        
    }

    return (EXIT_SUCCESS);
}





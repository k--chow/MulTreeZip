//*****************************************************************/
/*
This is Tree*Zip, a compression software for phylogenetic trees. 
It is based on HashBase, which is in turn based off of HashRF and HashCS,
developed by SeungJin Sul

(c) 2009 Tree*Zip: Suzanne Matthews
(c) 2009 HashBase: Suzanne Matthews
(c) 2009 HashRF : SeungJin Sul 
(c) 2009 HashCS : SeungJin Sul

  ***** Updated 2013 : Ruchi Chaudhary (UFL) for Mul-tree compression ****

*/
/*****************************************************/

#ifndef PARSING_H_
#define PARSING_H_

#include "global.h"
#include "hash.h"
#include "tree.h"
#include "tree_IO.h"
#include "tree_subtree_info.h"


#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

#ifndef NOHASH
#include <boost/unordered_map.hpp>
#endif

typedef boost::unordered_map<std::string, unsigned int> mul_cluster;

std::string decimal2binary(int number)
{
    string result = "";
    do
    {
        if ((number & 1)==0)
            result += "0";
        else
            result += "1";

        number >>= 1;
    } while(number);

    reverse(result.begin(), result.end());
    return result;
}

void initialize_hashtable(HashRFMap & vvec_hashrf, unsigned long long &M1, unsigned long long &M2) {    
    vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C);   
    M1 = vvec_hashrf._HF.getM1();
    M2 = vvec_hashrf._HF.getM2();
    vvec_hashrf._hashtab.resize(M1);    
}

    
boost::dynamic_bitset<> * dfs_compute_hash_w_bitstring(unsigned int node, aw::Tree &g_tree, aw::SubtreeParent<aw::Tree> &g_parent, aw::idx2name &g_taxa, HashRFMap &vvec_hashrf, int treeIdx, boost::dynamic_bitset<> * vect_bs[]) {
    if(g_tree.is_leaf(node)) {
        boost::dynamic_bitset<> * bs = new boost::dynamic_bitset<> (BITSET);
        std::string leaf_label = g_taxa[node];        
        int loc = guide[leaf_label];        
        (*bs)[loc] = 1;
        
        g_tree.set_hv1(node,vvec_hashrf._HF.getA1(lab2order[leaf_label]));
        g_tree.set_hv2(node,vvec_hashrf._HF.getA2(lab2order[leaf_label]));        
        vect_bs[node] = bs;        

        return bs;
    }
    else {
        std::vector<boost::dynamic_bitset<> *> ebs;
        int child = 0;        

        BOOST_FOREACH(const unsigned int &c,g_tree.children(node,g_parent.parent(node))) {
            child += 1;
            ebs.push_back(dfs_compute_hash_w_bitstring(c, g_tree, g_parent, g_taxa, vvec_hashrf, treeIdx, vect_bs));
        }

        // At this point, we find the bitset of this internal node
        boost::dynamic_bitset<> * bs = new boost::dynamic_bitset<> (BITSET);
        
        for (int i=0; i<child; ++i) {          
            if (ebs[i]) {
              bool carry = false, sum = false;
              for(int j=BITSET-1,jKK=0; j>=jKK; j--) {
                  if(!carry){
                      if((*bs)[j] && (*ebs[i])[j]) {
                        sum = false; carry = true;
                      }
                      else if((*bs)[j] || (*ebs[i])[j]) {
                        sum = true; carry = false;
                      }
                      else {
                          sum = false; carry = false;
                      }
                  }
                  else {
                      if((*bs)[j] && (*ebs[i])[j]) {
                        sum = true; carry = true;
                      }
                      else if((*bs)[j] || (*ebs[i])[j]) {
                        sum = false; carry = true;
                      }
                      else {
                          sum = true; carry = false;
                      }
                  }
                  (*bs)[j] = sum;
               }              
            }
            else {
              ERROR_exit("ERROR: null bitstring from leaf");              
            }
        }              
        vect_bs[node] = bs;

        unsigned long long temp1 = 0;
        unsigned long long temp2 = 0;

        unsigned long long m1 = vvec_hashrf._HF.getM1();
        unsigned long long m2 = vvec_hashrf._HF.getM2();

        BOOST_FOREACH(const unsigned int &c,g_tree.children(node,g_parent.parent(node))) {
          temp1 += g_tree.return_hv1(c);
          temp2 += g_tree.return_hv2(c);
        }

        unsigned long long temp;        
        temp = temp2 % m2;
        g_tree.set_hv2(node,temp);
        temp = temp1 % m1;
        g_tree.set_hv1(node,temp);

        BOOST_FOREACH(const unsigned int &c,g_tree.children(node,g_parent.parent(node))) {            
            g_tree.set_par_hv1(c,g_tree.return_hv1(node));
            g_tree.set_par_hv2(c,g_tree.return_hv2(node));
        }

        if(g_tree.root==node) {
            g_tree.set_par_hv1(node,NONODE);
            g_tree.set_par_hv2(node,NONODE);
        }

        return bs;
    }
    
}





#endif

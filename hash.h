//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		
//		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		HashRFMap: Class for hashing bitstrings

// ***** Updated 2013 : Ruchi Chaudhary (UFL) for Mul-tree compression ****
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).
//
//*****************************************************************/

#ifndef HASHRFMYHASH_HH
#define HASHRFMYHASH_HH

#include <map>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "hashfunc.h"
#include <bitset>

#ifndef NOHASH
#include <boost/unordered_map.hpp>
#endif

using namespace std;

typedef struct {
  float weights;
  unsigned long	long    _par_hv1;
  unsigned long	long    _par_hv2;
  unsigned long	long    _par_ord;  
} WG_PAR_STRUCT;
	
typedef struct {
  unsigned long	long    _hv2;
  boost::dynamic_bitset<> *_bs;
  bool is_leaf;
  bool once_par;   //true if at least one tree has this leaf as root's child
  vector<bool> par_root;   //true in a tree if this leaf is root's child
  vector<unsigned int> 	_vec_treeidx;
  vector<WG_PAR_STRUCT> wg_par;
  boost::unordered_map<unsigned int,unsigned int> _new_tree;  //store the last location of a tree in wg_par_pair  pair<treeIdx,lastLoc in wg_par_pair>  
} TREEIDX_STRUCT_T;

typedef vector<TREEIDX_STRUCT_T> V_BUCKET_T;
typedef vector<V_BUCKET_T> HASHTAB_T;
	

class HashRFMap {
    
    public:
      CHashFunc 	_HF;
      HASHTAB_T 	_hashtab;      // hash table that will store all hash values

      void uhashfunc_init(unsigned int t, unsigned int n, unsigned int c) {
          _HF.UHashfunc_init(t, n, c);
      }

      void hashrfmap_clear() {
          _hashtab.clear();
      }

      HashRFMap() {}
      ~HashRFMap() { }

      unsigned long long HashRFMap::hashing_bs(
                unsigned int treeIdx,
                bool leaf,
                bool par_r,
                unsigned long long hv1,
                unsigned long long hv2,
                unsigned long long par_hv1,
                unsigned long long par_hv2,
                unsigned long long par_ord,
                long double weight,
                boost::dynamic_bitset<> * tempbs)
        {
          
          unsigned sizeVec = _hashtab[hv1].size();

          if(!leaf || (leaf && par_r)) {
          std::cout<<"\n\n---------------------------";
          std::cout<<"\nTIdx:"<<treeIdx;
          std::cout<<"\nleaf:"<<leaf <<" par_r:"<<par_r;
          std::cout<<"\nhv1 & hv2:"<<hv1<<" "<<hv2;
          std::cout<<"\npar_hv1 & hv2:"<<par_hv1<<" "<<par_hv2;
          std::cout<<"\npar_ord:"<<par_ord;
          std::cout<<"\n"<<*tempbs;
          }

          if (sizeVec > 0) { //if the hash cell is not emply            
            for (unsigned int i=0; i<sizeVec; ++i) {
              if (_hashtab[hv1][i]._hv2 == hv2) {      //entering in a old bucket  - either cluster of a old tree or new tree

                if(*_hashtab[hv1][i]._bs != *tempbs)
                    ERROR_exit("DOUBLE Collision: start again...");

                WG_PAR_STRUCT wp;
                wp.weights = weight;
                wp._par_hv1 = par_hv1;
                wp._par_hv2 = par_hv2;
                wp._par_ord = par_ord;                
                _hashtab[hv1][i].wg_par.push_back(wp);
                unsigned int ordr = 0;
                
                              
                if (_hashtab[hv1][i]._vec_treeidx.back() != treeIdx) { //cluster of a new tree then update
                    unsigned int last_tree = _hashtab[hv1][i]._vec_treeidx.back();
                    _hashtab[hv1][i]._new_tree[treeIdx] = _hashtab[hv1][i]._new_tree[last_tree] + 1;
                    _hashtab[hv1][i]._vec_treeidx.push_back(treeIdx);
                    _hashtab[hv1][i].par_root.push_back(par_r);
                    if(par_r)
                        _hashtab[hv1][i].once_par=par_r;
                }
                else {  //same tree but another node ---mul-trees specific
                    if(par_r) {
                        _hashtab[hv1][i].once_par=par_r;
                        _hashtab[hv1][i].par_root.pop_back();
                        _hashtab[hv1][i].par_root.push_back(par_r);
                    }
                    _hashtab[hv1][i]._new_tree[treeIdx] = _hashtab[hv1][i]._new_tree[treeIdx] + 1;
                    unsigned int tot_trees = _hashtab[hv1][i]._vec_treeidx.size();
                    if(tot_trees==1)
                        ordr = _hashtab[hv1][i]._new_tree[treeIdx];
                    else {                        
                        unsigned int prev_tree = _hashtab[hv1][i]._vec_treeidx[tot_trees-2];
                        ordr = _hashtab[hv1][i]._new_tree[treeIdx] - _hashtab[hv1][i]._new_tree[prev_tree] - 1;
                    }
                }              
                                
                return ordr;
              }
            }
              //new bucket in the same hash cell
              TREEIDX_STRUCT_T bucket;
              bucket._hv2 = hv2;
              bucket._bs = tempbs;
              bucket.is_leaf = leaf;              
              bucket._vec_treeidx.push_back(treeIdx);
              bucket.once_par = false;
              bucket.par_root.push_back(par_r);
              if(par_r)
                  bucket.once_par=par_r;

              WG_PAR_STRUCT wp;
              wp.weights = weight;
              wp._par_hv1 = par_hv1;
              wp._par_hv2 = par_hv2;
              wp._par_ord = par_ord;              
              bucket.wg_par.push_back(wp);

              bucket._new_tree[treeIdx] = 0; //location of just entered element
              _hashtab[hv1].push_back(bucket);              

            return 0;
            
          }
          else if (sizeVec == 0) {//if there is a bipartition, but sizeVec is 0
            TREEIDX_STRUCT_T bucket; //create a new bucket, and add it to it
            bucket._hv2 = hv2;
            bucket._bs = tempbs;
            bucket.is_leaf = leaf;            
            bucket._vec_treeidx.push_back(treeIdx);
            bucket.once_par = false;
            bucket.par_root.push_back(par_r);
            if(par_r)
                bucket.once_par=par_r;
            //std::cout<<"\n--"<<bucket.once_par;

            WG_PAR_STRUCT wp;
            wp.weights = weight;
            wp._par_hv1 = par_hv1;
            wp._par_hv2 = par_hv2;
            wp._par_ord = par_ord;            
            bucket.wg_par.push_back(wp);            

            bucket._new_tree[treeIdx] = 0; //location of last cluster for treeIdx
            
            _hashtab[hv1].push_back(bucket);                            

            return  0;
         }

      }
};


#endif


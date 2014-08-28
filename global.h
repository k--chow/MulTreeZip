//*****************************************************************/
/*
This is Tree*Zip, a compression software for phylogenetic trees. 
It is based on HashBase, which is in turn based off of HashRF and HashCS,
developed by SeungJin Sul

(c) 2009 Tree*Zip: Suzanne Matthews
(c) 2009 FastHashRF: Suzanne Matthews
(c) 2009 HashRF : SeungJin Sul 
(c) 2009 HashCS : SeungJin Sul
 ***** Updated 2013 : Ruchi Chaudhary (UFL) for Mul-tree compression ****
*/
/*****************************************************/
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include "hash.h"
#include <boost/random.hpp>

#ifndef NOHASH
#include <boost/unordered_map.hpp>
#endif

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

unsigned int C			  = 200;
unsigned NUM_TREES                 = 0; // number of trees
unsigned int UNQ_TREES          = 0;
unsigned NUM_TAXA                	= 0; // number of taxa
unsigned int NUM_BIPART             = 0;
unsigned int BITSET                      = 0;    //numbr of taxa and their multiplicities
const unsigned long long NONODE = ULONG_MAX;
const unsigned int NOUINT = UINT_MAX;
bool WEIGHTED = true;
unsigned long long M1=0;
unsigned long long M2=0;

std::vector<unsigned int> unq_tr;   //shows that unq_num[i] = j means, jth tree in g_trees is ith unique
std::vector<unsigned int> sim2unq;   //sim2unq[i] = j, means ith tree g_trees is identical to jth tree in unq_num

unsigned long long BIG_M1;
int saved_char = 0;


typedef boost::unordered_map<std::string, unsigned int> bitstr_guide; // {name, right most bit}
bitstr_guide guide;

std::vector<std::string> order;    // {labels}   standard order

boost::unordered_map<unsigned int, unsigned int> guide2order;   //  each position points to the order position of the taxon it belongs to

std::vector<bool> multi_labeled;   // true if mul-tree, else false
bool all_singly = true;
bool all_mul = true;

unsigned int tot_sing = 0;  //# of singly labeled trees
unsigned int tot_mul  = 0;   //# of mul-trees

typedef boost::unordered_map<std::string, unsigned int> label2order;   // {name, order}
label2order lab2order;


#ifndef AW_RANDOMGEN
boost::mt19937 rng;
#endif

//extern unsigned int BITSETSZ = 0;  //bit for each GID with their multiplicities
#endif

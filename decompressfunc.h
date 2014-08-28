//*****************************************************************/
/*
/* MulTreeZip decompression *******
 * Ruchi Chaudhary (UFL) **********
*/
/*****************************************************/

#include <iostream>
#include <string>
#include "hash.h"
#include "tree_name_map.h"

#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>


using namespace std;

#ifndef _DECOMPRESSFUNC_H_
#define _DECOMPRESSFUNC_H_

//decode the bipartition
void decBipartition(const string &bipar, string &d_bipart) {
    char last;
    string num;
    bool n_flag = false;

    //std::cout<<"\nBipartition:::"<<bipar;
    string bipar1 = bipar + 'm';

    //processing bipartition
    for (unsigned i=0,iEE=bipar1.length(); i<iEE; ++i) {
        char ch = bipar1.at(i);

        if((last=='l' || last=='k') && !n_flag && !(ch<=57 && ch>=48))  //we expect a number after l & k
            ERROR_exit("ERROR!");

        if(last=='c' && !(ch<=57 && ch>=48) && !n_flag) { //settle for last 'ab' if not followed by number
            d_bipart += "01";
            last = ' ';
        }

        if(last=='a' && !(ch<=57 && ch>=48) && !n_flag && ch!='b') {  //settle for last 'a' if not followed by 'b'
            d_bipart += "0";
            last = ' ';
        }

        if(n_flag && (!(ch<=57 && ch>=48) || ch=='m')) {  //settle for non_number after number
            unsigned int cp = atoi(num.c_str());
            if(last =='l') {
                for(int j=0; j<cp; ++j)
                    d_bipart += '0';
            } else if(last == 'k') {
                for(int j=0; j<cp; ++j)
                    d_bipart += '1';
            } else if(last == 'c') {
                for(int j=0; j<cp; ++j)
                    d_bipart += "01";
            } else {                
                ERROR_exit("ERROR");
            }

            n_flag = false;
            num = ' ';
            last = ' ';
        }

        if(ch=='m')  break; 

        if(ch == 'l') {
            last = 'l';
        } else if(ch == 'k') {
            last = 'k';
        } else if(ch == 'a') {
            last = 'a';
        } else if(ch == 'b') {
            if(last=='a')
                last = 'c';
            else
                d_bipart += '1';
        } else {   //number
            num += ch;
            n_flag = true;
        }
    }


    if(d_bipart.size()!=BITSET) {
        std::cout<<"\nbipart: "<<bipar;
        std::cout<<"\n\nd_bipart"<<d_bipart;
        ERROR_exit("decoded bipartion size ("<<d_bipart.size()<<") is not equal to bitset size ("<<BITSET<<")!");
    }
    

}

//decoding tree numbers
//B1 = 11, C = 2
unsigned int compute(string tr) {    
    if(tr.size()==1)
        return tr.at(0) - 65;
    unsigned int f = tr.at(0) - 65;
    ostringstream s;
    s << f;
    string ss = s.str();
    ss = ss + tr.substr(1);
    return atoi(ss.c_str());
}

//decode mul-trees, parent info ....
//B5a6b.1C#   return - <15,a6b.1> <17,#>
void decodeMulTrees(const string &stp, std::vector<std::pair<unsigned int, string> > &ind) {
    //std::cout<<"\nstp: "<<stp;
    string tr;
    unsigned int tt;
    string par; //parent strip
    unsigned int pos = 0;
    char flag;
    
    for (unsigned i=0,iEE=stp.length(); i<iEE; ++i)  {
        char ch = stp.at(i);        
        if(ch<=90 && ch>=65) {  //CAPITAL LETTER            
            tr = ch;
            flag = 'C';
            if(!par.empty()) {
                ind.push_back(std::pair<unsigned int, string>(tt,par));
                par.clear();
            }
        } else if(ch<=57 && ch>=48 || ch == '.'){  //number or dot
            if(ch=='.')
                flag = 'T';
            if(flag=='C')
                tr += ch;
            if(flag=='S' || flag=='T')
                par += ch;
        }
        else if(ch<=122 && ch>=97 || ch=='#') {    //small letter
            if(!tr.empty()) {
                unsigned int num = compute(tr);
                if(pos!=0) {  //not first
                    tt = num + ind.back().first;                    
                }
                else {    //first
                    tt = num;     
                    pos++;
                }
                tr.clear();
            }
            if(!par.empty() && flag != 'T') {
                ind.push_back(std::pair<unsigned int, string>(tt,par));
                par.clear();
            }

            par += ch;
            flag = 'S';
        } else
            ERROR_exit("Wrong tree index!");
    }
    if(!par.empty())
        ind.push_back(std::pair<unsigned int, string>(tt,par));

    
//    std::cout<<" \nTrees-par: ";
//    for (unsigned i=0,iEE=ind.size(); i<iEE; ++i)
//       std::cout<<" <"<<ind[i].first<<","<<ind[i].second<<"> ";
}

//decode mul-trees, parent info ....
//b6 or b4.b#   return - 16 or 14.1
void process(const string &stp, string &ind) {
    string par;
    bool dot = false;    
    
    for (unsigned i=0,iEE=stp.length(); i<iEE; ++i)  {
        char ch = stp.at(i);
        if(ch<=122 && ch>=97 || ch == '#') {  //letter...
//            if(!dot && !par.empty()) {
//                ind = par;
//                par.clear();
//                dot = false;
//            }
            if(ch=='#')
               ind = "R";
            else {
                unsigned int ii = ch-97;
                ostringstream s;
                s << ii;
                par += s.str();
            }
            
        } else if(ch<=57 && ch>=48 || ch=='.'){     //number...
            par += ch;
            if(ch=='.')
                dot = true;
        }
        else
            ERROR_exit("Wrong tree index!");
    }
    if(!par.empty()) 
        ind = par;

    //std::cout<<" parents: ";    
    //std::cout<<"]";
}




//decode singly-labeled trees....
void decodeTrees(const string &stp, std::vector<unsigned int> &ind) {
    string tr;
    unsigned int idx = NOUINT;
    if(stp.length()==0)
        return;
    for (unsigned i=0,iEE=stp.length(); i<iEE; ++i)  {
        char ch = stp.at(i);
        if(ch<=90 && ch>=65) {
            if(!tr.empty()) {
                unsigned int num = compute(tr);
                if(idx!=NOUINT) {
                    ind.push_back(ind.back()+ num);
                    idx++;
                }
                else {
                    idx = 0;
                    ind.push_back(num);
                }
                tr.clear();
            }
            tr = ch;
        } else if(ch<=57 && ch>=48){
            tr += ch;            
        }
        else
            ERROR_exit("Wrong tree index!"<<ch<<" "<<stp);
    }
    if(!tr.empty()) {               
        unsigned int num = compute(tr);
        if(idx!=NOUINT) {
            ind.push_back(ind.back() + num);
            idx++;
        }
        else {
            idx = 0;
            ind.push_back(num);
        }        
    }

//    std::cout<<" Trees: ";
//    BOOST_FOREACH(const unsigned int id, ind)
//        std::cout<<id<<" ";
}

bool sort_bitstrings(const std::pair<boost::dynamic_bitset<>, unsigned> & a, const std::pair<boost::dynamic_bitset<>, unsigned int> &b) {
  return a.second > b.second;
}

//building the singly_labeled tree ............................................................................................................................................
void buildSTree(boost::dynamic_bitset<> &bs, std::vector<std::pair<boost::dynamic_bitset<>,unsigned int> > &my_bips,aw::Tree &tree,aw::idx2name &names) {

    boost::unordered_map<std::string, unsigned int> lab2node;
    unsigned int leaf_par[BITSET];  //leaf order is their order in the bs

    unsigned int rt = tree.new_node();
    tree.root = rt;    

    for(unsigned int i=0; i<BITSET; ++i) {
        if(bs[BITSET-1-i]==1) { //add a child leaf           
            string label = order[guide2order[i]];            
            unsigned int leaf_node = tree.new_node();
            tree.add_edge(rt,leaf_node);
            leaf_par[i] = rt;
            lab2node[label] = leaf_node;            
            names.insert(aw::idx2name::value_type(leaf_node, label));            
        }
    }    
    sort(my_bips.begin(), my_bips.end(), sort_bitstrings);
    for(unsigned int i=0,iEE=my_bips.size(); i<iEE; ++i) {
        if(my_bips[i].second == 1) {  //if leaf then child of root as encoded...
            continue;
        }
        boost::dynamic_bitset<> next_bs = my_bips[i].first;
        unsigned int new_node = tree.new_node();       
        unsigned int par = NONODE;
        
        for(int j=0; j<BITSET; ++j) {
            if(next_bs[BITSET-1-j]==1) {
                string label = order[guide2order[j]];
                unsigned int leaf_node = lab2node[label];               
                
                if(par==NONODE) {
                    par = leaf_par[j];                   
                    tree.add_edge(par,new_node);                  
                } else{                    
                    if(par != leaf_par[j])
                        ERROR_exit("ERROR!");
                }

                tree.remove_edge(par,leaf_node);                
                tree.add_edge(new_node,leaf_node);               
                leaf_par[j] = new_node;                
            }
        }
    }    
}


//building the multi_labeled tree ............................................................................................................................................
void buildMTree(std::vector<std::vector<string> > &ind, std::vector<std::pair<boost::dynamic_bitset<>,unsigned int> > &my_bips, aw::Tree &tree,aw::idx2name &names) {
    unsigned int rt = tree.new_node();
    tree.root = rt;

    boost::unordered_multimap<string, string> par2chd;    
    for(int i=0, iEE = my_bips.size(); i<iEE; ++i) {        
        if(my_bips[i].second != NOUINT ) {
            unsigned int pos_vec = my_bips[i].second;  //position in ind vector
            for(int j=0,jEE=ind[pos_vec].size(); j<jEE; ++j) {
                string par = ind[pos_vec][j];                
                if(par!="R" && par.find_first_of(".")==-1)
                    par = par + ".0";
                ostringstream ss;
                ss << i;
                string chld = ss.str();
                ss.str(""); ss << j;
                chld = chld + "." + ss.str();   // i.j
                par2chd.insert(std::pair<string, string>(par,chld));                
            }
        }
    }

    std::queue<std::pair<string,unsigned int> > que;   //string (reperesnting i.j i from my_bips and j from ind) and corresponding node in tree
    que.push(std::pair<string,unsigned int>("R",0));
    while(!que.empty()) {
        std::pair<string,unsigned int> elm = que.front();
        que.pop();
        string str = elm.first;
        unsigned int node = elm.second;

        {
            unsigned int dot = str.find_first_of(".");
            string ns = str.substr(0,dot);
            int gd1 = atoi(ns.c_str());
            boost::dynamic_bitset<> bp = my_bips[gd1].first;            
        }

        std::pair<boost::unordered_multimap<string, string>::iterator,boost::unordered_multimap<string, string>::iterator> ret = par2chd.equal_range(str);
        if(par2chd.find(str) == par2chd.end()) {    //If "NO CHILD" then add the leaf nodes under it and CONTINUE...
            int gd;
            unsigned int dot = str.find_first_of(".");
            string ns = str.substr(0,dot);
            gd = atoi(ns.c_str());
            boost::dynamic_bitset<> bp = my_bips[gd].first;            
            if(bp.count()>0) {
                //process it to attach all leaves
                for(unsigned int v=0,vEE=BITSET; v<vEE; ++v) {
                    if(bp[BITSET-1-v]==1) { //add a child leaf
                        string label = order[guide2order[v]];
                        unsigned int posi = guide[label];
                        int b_cps = posi - v;
                        int cnt = pow(2.0,b_cps);                        
                        for(int ct=0,ctEE=cnt; ct<ctEE; ct++) {                            
                            unsigned int leaf_node = tree.new_node();
                            tree.add_edge(node,leaf_node);
                            names.insert(aw::idx2name::value_type(leaf_node, label));
                        }
                    }
                }
            }
            else {
                ERROR_exit("Error!");
            }
            continue;
        }
        boost::dynamic_bitset<> bs(BITSET); //empty bit-string
        for (boost::unordered_multimap<string, string>::iterator it=ret.first; it!=ret.second; ++it) { //loop for each child 
            //computing sum of each child bitset with bs
            int gd;
            unsigned int dot = it->second.find_first_of(".");
            string ns = it->second.substr(0,dot);
            gd = atoi(ns.c_str());
            boost::dynamic_bitset<> bp = my_bips[gd].first;            
            if(node != tree.root) {
                bool carry = false, sum = false;
                for(int j=0,jKK=BITSET; j<jKK; j++) {
                  if(!carry){
                      if((bs)[j] && (bp)[j]) {
                        sum = false; carry = true;
                      }
                      else if(bs[j] || bp[j]) {
                        sum = true; carry = false;
                      }
                      else {
                          sum = false; carry = false;
                      }
                  }
                  else {
                      if(bs[j] && bp[j]) {
                        sum = true; carry = true;
                      }
                      else if(bs[j] || bp[j]) {
                        sum = false; carry = true;
                      }
                      else {
                          sum = true; carry = false;
                      }
                  }
                  bs[j] = sum;
                }
            }
            bool lf = false;
            if(bp.count()==1) {
                unsigned int fof = BITSET-1 - bp.find_first();
                string tx = order[guide2order[fof]];                
                if(fof==guide[tx])
                    lf = true;
            }            

            if(!lf) { //if not leaf node then process later
                const string &ptr = it->second;
                unsigned int nd = tree.new_node();
                tree.add_edge(node,nd);
                que.push(std::pair<string,unsigned int>(ptr,nd));
            } else {  //if it's leaf node then attach label etc
                unsigned int leaf_node = tree.new_node();
                tree.add_edge(node,leaf_node);
                unsigned int pos_one = BITSET-1 - bp.find_first();
                string label = order[guide2order[pos_one]];
                names.insert(aw::idx2name::value_type(leaf_node, label));
            }
        }        
        if(node != tree.root) {
            //subtract bs from bp for attaching additional leaves under node
            int gd;
            unsigned int dot = str.find_first_of(".");
            string ns = str.substr(0,dot);
            gd = atoi(ns.c_str());
            boost::dynamic_bitset<> bp = my_bips[gd].first;
            boost::dynamic_bitset<> bk(BITSET);

            {   //method for subtracting bitsets...
                
                bool borrow = false, sub = false;
                for(int j=0,jKK=BITSET; j<jKK; j++) {
                   
                  if(!borrow){
                      if(bp[j] && !bs[j]) {
                        sub = true; borrow = false;
                      }
                      else if(!bp[j] && bs[j]) {
                        sub = true; borrow = true;
                      }
                      else {
                          sub = false; borrow = false;
                      }
                  }
                  else {
                      if(!bp[j] && bs[j]) {
                        sub = false; borrow = true;
                      }
                      else if(bp[j] && !bs[j]) {
                        sub = false; borrow = false;
                      }
                      else {
                          sub = true; borrow = true;
                      }
                  }
                  bk[j] = sub;                  
                }
            }            

            if(bk.count()>0) {
                //process it to attach all leaves
                for(unsigned int v=0,vEE=BITSET; v<vEE; ++v) {
                    if(bk[BITSET-1-v]==1) { //add a child leaf
                        string label = order[guide2order[v]];
                        unsigned int posi = guide[label];
                        int b_cps = posi - v;
                        int cnt = pow(2.0,b_cps);                        
                        for(int ct=0,ctEE=cnt; ct<ctEE; ct++) {                            
                            unsigned int leaf_node = tree.new_node();
                            tree.add_edge(node,leaf_node);
                            names.insert(aw::idx2name::value_type(leaf_node, label));
                        }
                    }
                }
            }
        }
    }
}


//Decompress the mtz file...
void decompression(string in_file, std::vector<aw::Tree> &trees, std::vector<aw::idx2name> &names) {
    int ntaxa = 0;
    
    ifstream fin(in_file.c_str(), ios::binary);
    if (!fin) ERROR_exit("cannot open file!\n");

    string line, line_type;

    {   //READING FIRST LINE...
        bool saferead = getline(fin, line);
        if(!saferead)
            ERROR_exit("ERROR!");
        int pos = line.find_first_of(" ");
        line_type = line.substr(0, pos);
        if (line_type != "TX")
            ERROR_exit("Error!");
       
        line = line.substr(pos+1);  //string after 'TX'
        int posCR = line.find_first_of("\r");   //check \r (carriage return) for the newline in Mac & \r\n in windows
        line = line.substr(0,posCR);

        std::vector<std::string> segments;
        boost::split(segments, line, boost::is_any_of(" "));

        unsigned int last_guide = NOUINT;

        unsigned int c_ct = 0;
        //std::cout<<"\ncopies: ";
        BOOST_FOREACH(const string &c,segments) {
            string taxa_seg = c;            
            int cap_pos = taxa_seg.find_first_of("^");
            if(cap_pos==-1) {  //no duplicates of this taxa
                //std::cout<<" 1";
                c_ct += 1;
                order.push_back(taxa_seg);
                lab2order[taxa_seg] = order.size()-1;
                if(last_guide==NOUINT) {
                    last_guide = 0;
                    guide[taxa_seg] = last_guide;                    
                } else {
                    guide[taxa_seg] = last_guide + 1;
                    last_guide += 1;
                }                
            } else {  //duplicates then count them too..
                string copies = taxa_seg.substr(cap_pos+1);
                //std::cout<<" "<<copies;
                unsigned int cp = atoi(copies.c_str());
                c_ct += cp;
                //std::string bin_string = decimal2binary(cp);
                //int len_bin = bin_string.length();  //number of needed bits in binary

                taxa_seg = taxa_seg.substr(0,cap_pos);

                order.push_back(taxa_seg);
                if(last_guide==NOUINT) {
                    guide[taxa_seg] = cp - 1;
                    last_guide = cp - 1;
                }
                else {
                    guide[taxa_seg] = last_guide + cp;
                    last_guide += cp;
                }
                lab2order[taxa_seg] = order.size() - 1;

                
            }
            ntaxa++;
        }


        
        NUM_TAXA = ntaxa;
        BITSET = last_guide+1;
        //std::cout<<"\nBITSET"<<BITSET<<" "<<c_ct;
        

        for(unsigned int i=0; i<BITSET; ++i)
            guide2order[i] = NOUINT;

        //std::cout<<"\n";
        for(int i=0,iEE=order.size(); i<iEE; ++i) {
            string str = order[i];            
            unsigned int gd = guide[str];
            
            guide2order[gd] =  i;
            //std::cout<<" "<<str<<"<>"<<gd<<"<>"<<guide2order[gd];
        }
        
        unsigned int last;
        for(int i=BITSET-1; i>=0; --i) {
            if(guide2order[i]!=NOUINT)
                last = guide2order[i];
            else
                guide2order[i] = last;
        }       
    }

    {   //READING SECOND LINE...
        bool saferead = getline(fin, line);
        if(!saferead)
            ERROR_exit("ERROR!");
        
        //NUMBER OF TREES.....
        int pos_col = line.find_first_of(":");
        line_type = line.substr(0, pos_col);
        if (line_type != "NTR")
            ERROR_exit("Error!");
        line = line.substr(pos_col+1);  //string after "NTR"
        int pos_spac = line.find_first_of(" ");
        string num_t = line.substr(0,pos_spac);
        NUM_TREES = atoi(num_t.c_str());
        MSG("Total trees: "<<NUM_TREES);
        //std::cout<<"# Trees: "<<NUM_TREES;

        //NUMBER OF UNIQUE TREES.....
        line = line.substr(pos_spac+1);  //line number of trees
        pos_col = line.find_first_of(":");
        line_type = line.substr(0, pos_col);
        if (line_type != "UTR")
            ERROR_exit("Error!");
        line = line.substr(pos_col+1);  //string after "UTR"
        pos_spac = line.find_first_of(" ");
        num_t = line.substr(0,pos_spac);
        UNQ_TREES = atoi(num_t.c_str());
        MSG("Unique trees: "<<UNQ_TREES);
        //std::cout<<"# Unique Trees: "<<UNQ_TREES;

        // WEIGHTED or UNWEIGHTED...
        line = line.substr(pos_spac+1);
        pos_spac = line.find_first_of(" ");
        line_type = line.substr(0, pos_spac);
        if(line_type == "WG") {
            MSG("Weighted trees...");
            WEIGHTED = true;
        } else if(line_type == "NWG") {
            MSG("Unweighted trees...");
            WEIGHTED = false;
        } else
            ERROR_exit("Neigher WG nor NWG!");
        //std::cout<<" Weighted:"<<WEIGHTED;

        // ALL-MULs or ALL_SINGLY...
        line = line.substr(pos_spac+1);
        pos_spac = line.find_first_of(" ");
        line_type = line.substr(0, pos_spac);
        if(line_type == "SG") {
            all_mul = false;
            MSG("All singly-labeled trees...");
        } else if(line_type == "ML") {
            all_singly = false;
            MSG("All multi-labeled trees...");
        } else if(line_type == "SM") {
            MSG("Singly- and multi-labeled trees...");
            all_mul = false;
            all_singly = false;
        }
        else
            ERROR_exit("ERROR!");
        
        //NUMBER OF BIPARTITIONS.....
        line = line.substr(pos_spac+1);  //line number of trees
        pos_col = line.find_first_of(":");
        line_type = line.substr(0, pos_col);
        if (line_type != "NBP")
            ERROR_exit("Error!");
        line = line.substr(pos_col+1);  //string after "NBP"
        pos_spac = line.find_first_of("\r");
        string num_bipart = line.substr(0,pos_spac);
        NUM_BIPART = atoi(num_bipart.c_str());
        //std::cout<<"# Bipart: "<<NUM_BIPART;

    }    
    
    //reading bipartitions,trees....
    //NONWEIGHTED--------------------------------------------------------------------------------------------
    if(!WEIGHTED) {
        
        std::vector<std::vector<unsigned int> > sing_trs;  //indices of sing-trees
        std::vector<bool> plus;        
        bool is_mul[UNQ_TREES];
        
        std::vector<boost::dynamic_bitset<> > biparts; //bipartitions
        
        //unsigned int bips[NUM_BIPART][UNQ_TREES];         //each row - next bipart, column tree (mul or sing) indices

        int **bips = new int*[NUM_BIPART];
        for(int i = 0; i < NUM_BIPART; ++i) {
            bips[i] = new int[UNQ_TREES];
        }
        
        std::vector<std::vector<string> > vec;  //for mul-trees
        unsigned int count = 0;
        
        for (int i = 0; i < UNQ_TREES; ++i)
            is_mul[i] = false;
        
        for(int i=0,iEE=NUM_BIPART; i<iEE; ++i)
            for(int j=0,jEE=UNQ_TREES; j<jEE; ++j)
                bips[i][j] = NOUINT;

        std::vector<std::vector<std::pair<unsigned int, string> > > mul_inf;       //each row of mul_inf corresponds to a bipartition, storing mul-tree and string (parent,order) pairs
             
        while(getline(fin, line) && line.find("DUPS:")==-1) {            
            string bipar;
            string d_bipart;
            unsigned int next_pos;                  

            for (unsigned i=0,iEE=line.length(); i<iEE; ++i)  {
                char ch = line.at(i);
                if((ch<=57 && ch>=48) || ch==97 || ch==98 || ch==107 || ch==108)
                    bipar = bipar+ch;
                else{
                    next_pos = i;
                    break;
                }
            }

            //decode bipartition            
            decBipartition(bipar,d_bipart);            
            biparts.push_back(boost::dynamic_bitset<>(d_bipart));          
            
            //decode trees...
            line = line.substr(next_pos);
            string last_str;
            last_str = "\r";            

            std::vector<std::pair<unsigned int, string> > ind;
            if(all_mul) { //decode singly-labeled trees
                next_pos = line.find_first_of(last_str);
                if(next_pos == -1)
                    ERROR_exit("Error!");
                string mul_tr = line.substr(0,next_pos);
                if(mul_tr.size()!=0)
                    decodeMulTrees(mul_tr,ind);

                mul_inf.push_back(ind);
            } else if(all_singly) { //decode multi-labeled trees
                bool ps;
                if(line.at(0)=='+') {
                    ps = true;                    
                }
                else if(line.at(0)=='-') {
                    ps = false;                    
                }
                else
                    ERROR_exit("ERROR!");

                line = line.substr(1,line.size()-1);
                next_pos = line.find_first_of(last_str);
                if(next_pos == -1)
                    ERROR_exit("Error!");
                string sing_tr = line.substr(0,next_pos);                
                std::vector<unsigned int> sgs;
                decodeTrees(sing_tr,sgs);
                sing_trs.push_back(sgs);
                plus.push_back(ps);
            } else { //decode mul-trees then sing-trees..                
                bool ps = true;
                next_pos = line.find_first_of("+");
                if(next_pos == -1) {
                    next_pos = line.find_first_of("-");
                    ps = false;
                } if(next_pos == -1)
                    ERROR_exit("ERROR!");
                string mul_tr = line.substr(0,next_pos);               
                
                if(mul_tr.size()!=0) {                    
                    decodeMulTrees(mul_tr,ind);
                }

                mul_inf.push_back(ind);
                line = line.substr(next_pos+1);
                next_pos = line.find_first_of(last_str);
                if(next_pos == -1)
                    ERROR_exit("Error!");
                string sing_tr = line.substr(0,next_pos);
                
                std::vector<unsigned int> sgs;
                decodeTrees(sing_tr,sgs);
                sing_trs.push_back(sgs);
                plus.push_back(ps);
                //std::cout<<" here22";
            }
           
            if(!all_singly) {              
                //process ind vector for each tree in this cluster
                unsigned int next_t = NOUINT;
                std::vector<string> pars;
                for(unsigned int i=0,iEE=ind.size(); i<iEE; ++i) {
                    unsigned int tr = ind[i].first;

                    if(tr>=UNQ_TREES)
                        ERROR_exit("Wrong tree number!"<<tr<<" "<<UNQ_TREES);

                    string strip = ind[i].second;                    

                    is_mul[tr] = true;
                    string par;
                    process(strip,par);
                    
                    if(next_t==NOUINT) {
                        next_t = tr;
                        pars.push_back(par);
                    }
                    else {
                        if(tr==next_t)
                            pars.push_back(par);
                        else {
                            vec.push_back(pars);
                            bips[count][next_t] = vec.size() - 1;
                            pars.clear();
                            pars.push_back(par);
                            next_t = tr;
                        }
                    }
                }
                if(pars.size()!=0) {
                    vec.push_back(pars);
                    bips[count][next_t] = vec.size() - 1;
                }                
            }
            count++;           
        }
             
        //process bips for singly-labeled trees...
        if(!all_mul) {            
            for(int i=0,iEE=NUM_BIPART; i<iEE; ++i) {
                std::vector<unsigned int> vec = sing_trs[i];
                bool fg = plus[i];
                if(fg) {
                    for(int j=0,jEE=vec.size(); j<jEE; ++j)
                        bips[i][vec[j]] = 1;
                }
                else {
                    for(int j=0,jEE=UNQ_TREES; j<jEE; ++j)
                        if(!is_mul[j])
                            bips[i][j] = 1;
                    for(int j=0,jEE=vec.size(); j<jEE; ++j)
                        bips[i][vec[j]] = NOUINT;
                }
            }
        }

        for(int i=0,iEE=UNQ_TREES; i<iEE; ++i) 
            multi_labeled.push_back(is_mul[i]);            
        
        //std::cout<<"\nBefore Building Trees!";
      
        std::vector<aw::Tree> u_trees;
        std::vector<aw::idx2name> u_names;
        
        //Building trees now ...
        for(int tr=0; tr<UNQ_TREES; ++tr) {   // runs number of uniq tree times            
            if(multi_labeled[tr]) {                
                //build ML tree...
                std::vector<std::pair<boost::dynamic_bitset<>,unsigned int> > my_bips;   //tr column of bips with bipartitions
                for(int row=0; row<NUM_BIPART; ++row) {
                     my_bips.push_back(std::pair<boost::dynamic_bitset<>,unsigned int>(biparts[row],bips[row][tr]));                     
                }

                aw::Tree tree;
                aw::idx2name name;
                buildMTree(vec,my_bips,tree,name);
                u_trees.push_back(tree);
                u_names.push_back(name);               
            }
            else {
                boost::dynamic_bitset<> bs(BITSET);
                for(int row=0; row<NUM_BIPART; ++row) {                   
                    if(bips[row][tr]!=NOUINT)
                        bs |= biparts[row];                    
                }
                
                std::vector<std::pair<boost::dynamic_bitset<>,unsigned int> > my_bips;
                for(int row=0; row<NUM_BIPART; ++row) {
                    if(bips[row][tr]!=NOUINT) {
                        my_bips.push_back(std::pair<boost::dynamic_bitset<>,unsigned int>(biparts[row],biparts[row].count()));
                    }
                }
                
                aw::Tree tree;
                aw::idx2name name;
                buildSTree(bs, my_bips, tree, name);
                u_trees.push_back(tree);
                u_names.push_back(name);
            }
        }

        unsigned int all_trees[NUM_TREES];  //pointer to unique trees
        bool unq_done[UNQ_TREES];
        for(int i=0; i<NUM_TREES; ++i)
            all_trees[i] = NOUINT;
        for(int i=0; i<UNQ_TREES; ++i)
            unq_done[i] = false;
        
        while(getline(fin, line)) {
            unsigned int ed = line.find_first_of("\r");
            line = line.substr(0,ed);
            
            std::vector<unsigned int> positions;
            decodeTrees(line,positions);
            unq_done[positions[0]] = true;
            for(int i=1,iEE=positions.size(); i<iEE; ++i)
                all_trees[positions[i]] = positions[0];
        }

        int last = 0;
        for(int i=0; i<UNQ_TREES; ++i) {
            if(!unq_done[i]) {
                for(int j=last; j<NUM_TREES; ++j) {
                    if(all_trees[j]==NOUINT) {
                        all_trees[j] = i;
                        last = j+1;
                        break;
                    }
                }
            }
        }

        //storing trees....
        for(int i=0; i<NUM_TREES; ++i) {
            unsigned int which_tr = all_trees[i];            
            trees.push_back(u_trees[which_tr]);
            names.push_back(u_names[which_tr]);
        }

        //clean up bips.....
        for(int i = 0; i < NUM_BIPART; ++i) {
            delete [] bips[i];
        }
        delete [] bips;
    }




}

#endif

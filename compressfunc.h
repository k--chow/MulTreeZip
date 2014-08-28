//*****************************************************************/
/*
 Encoding hash file : by Ruchi Chaudhary

*/
/*****************************************************/

#include <iostream>
#include <string> 
#include "hash.h"
#include "tree_name_map.h"

#include <boost/algorithm/string.hpp>


using namespace std;

#ifndef _COMPRESSFUNC_H_
#define _COMPRESSFUNC_H_

//function for efficient encoding weight
//tree wg = 234 then code = x4
string codeWeight(const float &wg) {
    char smallAlpha[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
    string final_string;

    ostringstream sd;
    sd << wg;
    string wg_str = sd.str();    
    int loc_dot = wg_str.find(".");
    string dec_str;
    if(loc_dot==-1)
        dec_str = "";
    else
        dec_str = wg_str.substr(loc_dot);
    
    int int_part = (int)wg;       
  
    if(int_part < 26){
        final_string = smallAlpha[int_part];
    }
    else {
        ostringstream s;
        s << int_part;
        string ss = s.str();        
        string s_sub = ss.substr(0,2);        
        unsigned int x = atoi(s_sub.c_str());        

        if(x<26){            
            final_string = smallAlpha[x] + ss.substr(2,ss.length());
        } else {            
            s_sub = ss[0];
            x = atoi(s_sub.c_str());
            final_string = smallAlpha[x] + ss.substr(1,ss.length());
        }
    }
    
    return final_string + dec_str;
}

//function for efficient encoding node parent
//parent par = 234 then code = X4
string codeParent(const unsigned int &par) {
    char smallAlpha[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
    string final_string;

    if(par < 26){
        final_string = smallAlpha[par];
    }
    else {
        ostringstream s;
        s << par;
        string ss = s.str();        
        string s_sub = ss.substr(0,2);        
        unsigned int x = atoi(s_sub.c_str());        

        if(x<26){            
            final_string = smallAlpha[x] + ss.substr(2,ss.length());
        } else {            
            s_sub = ss[0];
            x = atoi(s_sub.c_str());
            final_string = smallAlpha[x] + ss.substr(1,ss.length());
        }
    }    
    return final_string;
}

//function for efficient encoding tree numbers
//tree nb = 234 then code = X4
string codeNumb(const unsigned int &nb) {
    char bigAlpha[] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};    
    string final_string;
    
    if(nb < 26){
        final_string = bigAlpha[nb];
    }
    else {        
        ostringstream s;
        s << nb;
        string ss = s.str();
        //cout<<" *ss "<<ss;
        string s_sub = ss.substr(0,2);
        //cout<<" *s_sub "<<s_sub;
        unsigned int x = atoi(s_sub.c_str());
        //cout<<" *x "<<x;

        if(x<26){
            //cout<<" *REST "<<ss.substr(2,ss.length());
            final_string = bigAlpha[x] + ss.substr(2,ss.length());
        } else {
            //cout<<" *REST "<<ss.substr(1,ss.length());
            s_sub = ss[0];
            x = atoi(s_sub.c_str());
            final_string = bigAlpha[x] + ss.substr(1,ss.length());
        }
    }
    //std::cout<<"\n"<<final_string;
    return final_string;
}

//formatting bitstring: Single 0 = a, 1 = b,  Multiple 0 = l, 1 = k
// ababab = ab3 &&  0000 = l4 && 1111 = k4
void runlength(boost::dynamic_bitset<> * temp, string &rlstring) {
    
    unsigned int curr = (*temp)[0];
    string tempstr;
    ostringstream s;
    unsigned int currcnt = 1;
    unsigned int newvalb = 0;
    unsigned int pairBA = 0;
    for (unsigned int j=1; j<(*temp).size(); ++j) {
        //cout<<j;
        newvalb = (*temp)[j];
        //cout<<" "<<newvalb;

        if (newvalb != curr) {
            if( currcnt > 1) {
                if (curr == 1) {
                    s.str("");
                    s << currcnt;
                    tempstr = tempstr + "k" + s.str();
                }
                else {
                    s.str("");
                    s << currcnt;
                    tempstr = tempstr + "l" + s.str();
                }
            } else {
                if (curr == 1)  tempstr = tempstr + "b";
                else  tempstr = tempstr + "a";
            }
            curr = newvalb;
            currcnt = 1;
       }
       else{  currcnt++;
       }
    }

    if (currcnt > 1) {
       if (curr == 1)  {
           s.str("");
           s << currcnt;
           tempstr = tempstr + "k" + s.str();
       }
       else {
           s.str("");
           s << currcnt;
           tempstr = tempstr + "l" + s.str();
       }
    } else{
       if(curr == 1)  tempstr = tempstr + "b";
       else  tempstr = tempstr + "a";
    }

    //cout<<"\n"<<tempstr;
    int vv = 0;
    string buffer;
    int tmpstrlen = tempstr.length();
    for (int i=0; i<tmpstrlen; ++i) {
        char ch = tempstr[i];        
        if(ch=='a') {
            if(vv==0) {
                vv = 1; buffer += 'a';
            }
            else {
                int buf_len = buffer.length();
                if((buf_len % 2) == 0) ERROR_exit("ERROR");
                if((buf_len-1)/2 > 1) {
                    int llen = (buf_len-1)/2;
                    ostringstream s;
                    s << llen;
                    rlstring += "ab" + s.str() + "a";
                }
                else
                    rlstring += "aba";
                buffer = "a";
            }
        } else if(ch=='b') {
            if(vv==1) {
                vv = 0;
                buffer += 'b';
            } else {
                int buf_len = buffer.length();
                if(buf_len!=0){
                    if((buf_len-1)/2 > 1){
                        ostringstream s;
                        s << (buf_len-1)/2;
                        rlstring += "ab" + s.str() + "b";
                    }
                    else
                        rlstring += "abb";
                    buffer = "";
                } else {
                    rlstring += ch;
                }
            }
        } else {
            vv = 0;
            int buf_len = buffer.length();
            if(buf_len!=0) {
                int mod = buf_len % 2;
                ostringstream s;
                if(mod == 0) {
                    if (buf_len > 2) {
                        s << buf_len/2;
                        rlstring += "ab" + s.str();
                    } else
                        rlstring += "ab";
                }
                else {
                    if (buf_len > 3) {
                        s << (buf_len-1)/2;
                        rlstring += "ab" + s.str() + "a";
                    } else {
                        rlstring += buffer;
                    }
                }
                buffer = "";
            }
            rlstring += ch;            
        }
        int buf_len = buffer.length();
        if(i == tmpstrlen-1 && buf_len>0) {            
            int mod = buf_len % 2;
            ostringstream s;
            if(mod == 0) {
                if (buf_len > 2) {
                    s << buf_len/2;
                    rlstring += "ab" + s.str();
                } else
                    rlstring += "ab";
            }
            else {
                if (buf_len > 3) {
                    s << (buf_len-1)/2;
                    rlstring += "ab" + s.str() + "a";
                } else {
                    rlstring += buffer;
                }
            }            
        }
    }
    if(tempstr!=rlstring)
        saved_char += tempstr.length() - rlstring.length();        
}
 
//compresses a tree file into compressed form
void compressf(string out_file, HashRFMap &vvec_hashrf) {

    std::ofstream ouput_fs;
    if (!out_file.empty()) {
        ouput_fs.open(out_file.c_str());
        if (!ouput_fs) ERROR_exit("cannot write file '" << out_file << "'");
    }
    std::ostream &output = out_file.empty() ? std::cout : ouput_fs;

    output<<"TX ";

    
    for(unsigned int i=0,iEE=order.size(); i<iEE; ++i) {
        std::string name = order[i];
        int copies;        
        if(i==0) {
            copies = guide[name] + 1;
        }
        else {
            output<<" ";
            std::string last_name = order[i-1];
            copies = guide[name] - guide[last_name];
        }
        
        if(copies>1)
            output<<name<<"^"<<copies;
        else
            output<<name;
        
    }
    
    output<<endl;
    
    output<<"NTR:"<<NUM_TREES;
    output<<" UTR:"<<unq_tr.size();

    if(WEIGHTED)
        output<<" WG";
    else
        output<<" NWG";

    if(all_singly)
        output<<" SG";
    else if(all_mul)
        output<<" ML";
    else
        output<<" SM";

    //count up the number of unique bipartitions and compute bitstring order
    vector<string>  ordered_bitstrings;
    boost::unordered_map<string,unsigned int> hv1_hv2;
    boost::unordered_map<string,unsigned int> bitstring_ordering;
    unsigned int bipart_counter = 0;
    
    {   //bipartitions ordering 
        ostringstream ss;
        string loc;
        unsigned int ordering = 0;
        for(unsigned int hti =0,htiEE=vvec_hashrf._hashtab.size(); hti<htiEE; ++hti) {
            unsigned int sizeVec = vvec_hashrf._hashtab[hti].size();
            if(sizeVec>0) {
                //bipart_counter += sizeVec;
                for (unsigned int i=0; i<sizeVec; ++i) {
                    bipart_counter++;
                    if(!WEIGHTED) {   //we don't need to add leaves when their parents are not root for nonweighted
                        bool is_leaf = vvec_hashrf._hashtab[hti][i].is_leaf;
                        bool par_r = vvec_hashrf._hashtab[hti][i].once_par;
                        //std::cout<<"\n"<<is_leaf<<" "<<par_r;
                        if(is_leaf && !par_r) {
                            bipart_counter--;
                            continue;
                        }
                    } 

                    unsigned long long hv2 = vvec_hashrf._hashtab[hti][i]._hv2;

                    ss.str("");  ss << hti;
                    loc = ss.str();
                    ss.str("");  ss << i;
                    loc = loc + "." + ss.str();
                    
                    ordered_bitstrings.push_back(loc);
                    bitstring_ordering[loc] = ordering;                    

                    ss.str(""); ss << hti;
                    loc = ss.str();
                    ss.str(""); ss << hv2;
                    loc = loc + "." + ss.str();                    

                    hv1_hv2[loc] = ordering++;

                    //TESTING...
//                    boost::dynamic_bitset<> * temp = vvec_hashrf._hashtab[hti][i]._bs;
//                    cout<<"\nBitstring ";
//                    for (unsigned int j=0; j<BITSET; ++j) {
//                        cout<<" "<<(*temp)[j];
//                    }
                } //for everything in sizevec
            }
        }
    }
    
    output<<" NBP:"<<bipart_counter<<"\n";

    //start writing bitstrings with details in the file
    unsigned int x = 0, y = 0, p = 0;
    string place, xval;    
    
    BOOST_FOREACH(const string &c,ordered_bitstrings) {        
        place = c;

        //fatch hash-cell('x') and bucket('y') for this bitstring
        p = place.find_first_of(".");
        xval = place.substr(0,p);
        place = place.substr(p+1);
        x = atoi(xval.c_str());
        y = atoi(place.c_str());

        //printing formatted bitstring
        string ourstring;        
        runlength(vvec_hashrf._hashtab[x][y]._bs, ourstring);       
        output<<ourstring;

        bool this_leaf = false;
        if((*vvec_hashrf._hashtab[x][y]._bs).count()==1)
            this_leaf = true;

        //write the trees & parents
        unsigned int vec_ptrs = vvec_hashrf._hashtab[x][y]._vec_treeidx.size();

        //extract muls and singly labeled trees
        int countSing = 0, countMul = 0;
        vector<unsigned int> muls;        //stores indices in vec_tree_idx that have mul-trees
        vector<unsigned int> no_muls;   //stores indices in vec_tree_idx that have singly-labeled trees
        if(all_singly) {
            countSing = vec_ptrs;
            for (unsigned int j = 0; j < vec_ptrs; ++j) {
                no_muls.push_back(j);
            }
        }
        else {
            for (unsigned int j = 0; j < vec_ptrs; ++j) {
                unsigned int vec_ptr = j;
                unsigned int tr = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr];
                if(multi_labeled[tr]){
                    countMul++;
                    muls.push_back(vec_ptr);
                }
                else {
                    countSing++;
                    no_muls.push_back(vec_ptr);
                }
            }
        }
          
        if(!all_singly) { //Write MUL-TREES first...
            unsigned int ft = NOUINT;  //denote what we wrote last
            for (unsigned int k=0; k<countMul; ++k) {
                unsigned int vec_ptr = muls[k];
                unsigned int treeIdx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr];
                //std::cout<<"\n&&& "<<vvec_hashrf._hashtab[x][y].is_leaf<<vvec_hashrf._hashtab[x][y].par_root[vec_ptr];
                if(vvec_hashrf._hashtab[x][y].is_leaf && !vvec_hashrf._hashtab[x][y].par_root[vec_ptr])
                    continue;

                //write runlength treeIdx & compute its copies
                string treeCode;
                unsigned int loc2, loc1 = 0;
                if(vec_ptr != 0) {
                    unsigned int prvIndx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr-1];
                    loc1 = vvec_hashrf._hashtab[x][y]._new_tree[prvIndx] + 1;
                }

                if(ft!=NOUINT) {
                    unsigned int prevmulTree = vvec_hashrf._hashtab[x][y]._vec_treeidx[muls[ft]];
                    treeCode = codeNumb(treeIdx - prevmulTree);
                }
                else {
                    treeCode = codeNumb(treeIdx);
                }
                loc2 = vvec_hashrf._hashtab[x][y]._new_tree[treeIdx];
                output<<treeCode;
                ft = k;
                
                //compute parent's row.column
                for (unsigned int clst=loc1; clst<=loc2; ++clst) {
                    unsigned long long p_hv1 = vvec_hashrf._hashtab[x][y].wg_par[clst]._par_hv1;
                    unsigned long long p_hv2 = vvec_hashrf._hashtab[x][y].wg_par[clst]._par_hv2;
                    unsigned long long p_ord = vvec_hashrf._hashtab[x][y].wg_par[clst]._par_ord;

                    if(p_hv1==BIG_M1) {
                        if(!vvec_hashrf._hashtab[x][y].is_leaf || (vvec_hashrf._hashtab[x][y].is_leaf && vvec_hashrf._hashtab[x][y].par_root[vec_ptr]))
                            output<<"#";
                        continue;
                    }
                    if(vvec_hashrf._hashtab[x][y].is_leaf)
                        continue;

                    ostringstream ss1;
                    string location;
                    ss1 << p_hv1;
                    location = ss1.str();
                    ss1.str("");
                    ss1 << p_hv2;
                    location = location + "." + ss1.str();
                    unsigned int parent_row = hv1_hv2[location];                    
                                    
                    output<<codeParent(parent_row);
                    if(p_ord>0)
                        output<<"."<<codeParent(p_ord);
                }
                
            }
        }  //End-of writing MUL-TREES...

        if(!all_mul) {    //write SINGLY-LABELED trees            
            if((!WEIGHTED && countSing >= (60.0*(NUM_TREES-countMul))/100.0) && !(!WEIGHTED && this_leaf)) {
                //write trees that are not present....
                unsigned int prv = 0;
                vector<unsigned int> not_in;
                for (unsigned int ii=0; ii<vec_ptrs; ++ii) {
                    unsigned int id = vvec_hashrf._hashtab[x][y]._vec_treeidx[ii];

                    if(ii==0) {
                        if(id>0) {
                            for (unsigned int jj=0; jj<id; ++jj)
                                not_in.push_back(jj);
                        }
                    }
                    else {
                        for (unsigned int jj=prv+1; jj<id; ++jj)
                            not_in.push_back(jj);
                    }
                    prv = id;
                }
                unsigned int lastid = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptrs-1];
                if(lastid<(NUM_TREES-1))
                    for (unsigned int jj=lastid+1; jj<NUM_TREES; ++jj)
                        not_in.push_back(jj);

                //actual write down...
                output<<"-";
                for (unsigned int k=0; k<not_in.size(); ++k) {
                    unsigned int treeIdx = not_in[k];
                    string treeCode;
                    if(k!=0) {
                        unsigned int prevTree = not_in[k-1];
                        treeCode = codeNumb(treeIdx - prevTree);
                    }
                    else {
                        treeCode = codeNumb(treeIdx);
                    }
                    output<<treeCode;
                }
            }
            else {
                //write trees that are present...
                output<<"+";
                unsigned int ft = NOUINT;  //denote what we wrote last
                for (unsigned int k=0; k<countSing; ++k) {
                    unsigned int vec_ptr = no_muls[k];
                    unsigned int treeIdx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr];
                    //std::cout<<" *"<<vvec_hashrf._hashtab[x][y].par_root[vec_ptr];
                    if(vvec_hashrf._hashtab[x][y].is_leaf && !vvec_hashrf._hashtab[x][y].par_root[vec_ptr])
                        continue;
                    
                    string treeCode;
                    if(ft!=NOUINT) {
                        unsigned int prevmulTree = vvec_hashrf._hashtab[x][y]._vec_treeidx[no_muls[ft]];
                        treeCode = codeNumb(treeIdx - prevmulTree);
                    }
                    else {
                        treeCode = codeNumb(treeIdx);
                    }
                    output<<treeCode;
                    ft = k;
                }
            }
        }   //End-of writing SINGLY-LABELED trees
           
        if(WEIGHTED) {      //Add weights...
            output<<" ";            
            if(!all_singly) { //Write weights for Mul-trees...
                for (unsigned int k=0; k<countMul; ++k) {
                    unsigned int vec_ptr = muls[k];
                    unsigned int treeIdx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr];
                    unsigned int loc1 = 0,loc2;
                    if(vec_ptr!=0) {
                        unsigned int prvIndx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr-1];
                        loc1 = vvec_hashrf._hashtab[x][y]._new_tree[prvIndx] + 1;
                    }

                    loc2 = vvec_hashrf._hashtab[x][y]._new_tree[treeIdx];
                    
                    for (unsigned int clst=loc1; clst<=loc2; ++clst) {                        
                        float wg = vvec_hashrf._hashtab[x][y].wg_par[clst].weights;
                        string wgcode = codeWeight(wg);
                        output<<wgcode;                        
                    }
               }
            }
            if(!all_mul) {  //write weights for Singly-labeled trees...
                for (unsigned int k=0; k<countSing; ++k) {
                    unsigned int vec_ptr = no_muls[k];
                    unsigned int loc2;
                    unsigned int treeIdx = vvec_hashrf._hashtab[x][y]._vec_treeidx[vec_ptr];

                    loc2 = vvec_hashrf._hashtab[x][y]._new_tree[treeIdx];
                    float wg = vvec_hashrf._hashtab[x][y].wg_par[loc2].weights;
                    string wgcode = codeWeight(wg);
                    output<<wgcode;                                     
                }
            }
        } //End-of-writing weights...
        
        output<<"\n";
    }

    if(!WEIGHTED) {
        boost::unordered_multimap<unsigned int, unsigned int> Q;   //pair <uniq_pos,g_trees_pos>
        unsigned int num = 0;
        bool flag = false;
        
        BOOST_FOREACH(const unsigned int &c, sim2unq) {            
            if(c!=num)
                flag = true;            
            Q.insert(std::pair<unsigned int, unsigned int>(c,num));            
            num++;
        }
        if(flag)
            output<<"DUPS:\n";

        for(int c=0,cEE=unq_tr.size(); c<cEE; ++c) {
            std::pair<boost::unordered_multimap<unsigned int, unsigned int>::iterator,boost::unordered_multimap<unsigned int, unsigned int>::iterator> ret = Q.equal_range(c);
            unsigned int ent = 0;
            string treeCode = codeNumb(c);
            unsigned int prev = c;            
            for(boost::unordered_multimap<unsigned int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {                
                ent++;
                unsigned int pos = it->second;
                treeCode = treeCode + codeNumb(pos-prev);
                prev = pos;
            }
            if(ent>1)
                output<<treeCode<<"\n";           
        }
    }
}


//void decompress(string file);

#endif


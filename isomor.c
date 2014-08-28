/*
 * Implementation of mul-trees isomorphism algorithm described in
 * Ganapathy Ganeshkumar et al. "Pattern Identification in Biogeography"
 * IEEE/ACM TCBB Vol 3, 2006.
 *
 * Copyright (C) 2013 Ruchi Chaudhary
 * 
 */

#ifndef ISO_C_
#define ISO_C_

#include <stdlib.h>
#include <list>
#include <queue>
#include "tree.h"
#include "tree_name_map.h"

#include <boost/foreach.hpp>

#ifndef NOHASH
#include <boost/unordered_set.hpp>
#endif

#ifndef NOHASH
#include <boost/unordered_map.hpp>
#endif

//Laxicographically sort a list of n tuples of length k
//ptrs    list
//0,2,4   1,2,4,3,5,1 [tuples - (1,2),(4,3),(5,1)]
void sort_laxiK(std::vector<unsigned int> &ptrs, std::vector<unsigned int> &list, std::vector<unsigned int> &s_ptrs) {

    std::queue<unsigned int> QUEUE; //pointers to tuples in list
    int k = ptrs[1] - ptrs[0];
    int n = ptrs.size();

    for(int i=0; i<n; ++i)
        QUEUE.push(ptrs[i]);

    for(int j = k-1; j >= 0; --j) {
        boost::unordered_multimap<unsigned int, unsigned int> Q; //pairs elm, pointers
        unsigned int max =0;
        while(QUEUE.size()>0) {
            unsigned int ptr = QUEUE.front();
            QUEUE.pop();
            unsigned int aij = list[ptr+j];
            if(aij>max) max = aij;

            Q.insert(std::pair<unsigned int,unsigned int>(aij,ptr));
        }
        for(int jj = 0; jj <= max; ++jj) {
            std::pair<boost::unordered_multimap<unsigned int, unsigned int>::iterator,boost::unordered_multimap<unsigned int, unsigned int>::iterator> ret = Q.equal_range(jj);
            for (boost::unordered_multimap<unsigned int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
                const unsigned int &ptr = it->second;
                QUEUE.push(ptr);
            }
        }
    }

    for(int i=0; i<n; ++i) {
        unsigned int pt = QUEUE.front();
        QUEUE.pop();
        s_ptrs.push_back(pt);
    }
}


//Laxicographically sort a list of n tuples of varing length
//ptrs    list
//0,4,7   3,1,2,3,2,4,3,1,1 [tuples - (1,2,3),(4,3),(1)]
void sort_laxi(std::vector<unsigned int> &ptrs, std::vector<unsigned int> &lists, std::vector<unsigned int> &s_ptrs) {
    int n = ptrs.size(); //n - number of tuples
    std::vector<unsigned int> pairs;
    std::vector<unsigned int> ps;

    for(int i=0; i<n; ++i) {
        unsigned int ptr = ptrs[i];
        unsigned int elms = lists[ptr];

        unsigned int m = 0;
        for(int i=ptr+1; i<=ptr+elms; ++i) {
            ps.push_back(pairs.size());
            pairs.push_back(m);
            pairs.push_back(lists[i]);
            m++;
        }
    }
    std::vector<unsigned int> s_pairs;
    sort_laxiK(ps,pairs,s_pairs);

    std::vector<std::list<unsigned int> > NONEMPTY; //pairs elm, pointers
    unsigned int th = 0;

    std::list<unsigned int> elms;
    BOOST_FOREACH(const unsigned int &c,s_pairs) {
        int pos = pairs[c];
        unsigned int val = pairs[c+1];
        if(pos==th)
            elms.push_back(val);
        else {
            elms.sort();
            elms.unique();
            NONEMPTY.push_back(elms);
            elms.clear();
            elms.push_back(val);
            th = pos;
        }
    }
    //last round
    elms.sort();
    elms.unique();
    NONEMPTY.push_back(elms);

    std::vector<int> l;    //length of tuples
    int maxL = 0;
    boost::unordered_multimap<int, unsigned int> LENGTH; //pairs elm, pointers
    BOOST_FOREACH(const unsigned int &c,ptrs) {
        l.push_back(lists[c]);
        if(lists[c]>maxL) maxL = lists[c];
        LENGTH.insert(std::pair<int, unsigned int>(lists[c],c));
    }

    std::deque<unsigned int> QUEUE; //pointers to tuples in list
    for(int i = maxL; i > 0; --i) {

        //putting length i tuples to the QUEUE
        std::pair<boost::unordered_multimap<int, unsigned int>::iterator,boost::unordered_multimap<int, unsigned int>::iterator> ret = LENGTH.equal_range(i);
        std::deque<unsigned int>::iterator itr = QUEUE.begin();

        for (boost::unordered_multimap<int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
            const unsigned int &ptr = it->second;
            itr = QUEUE.insert(itr,ptr);
            itr++;
        }

        boost::unordered_multimap<unsigned int, unsigned int> Q; //pairs elm, pointers

        while(QUEUE.size()>0) {
            unsigned int ptr = QUEUE.front();
            QUEUE.pop_front();
            unsigned int aij = lists[ptr+i];
            Q.insert(std::pair<unsigned int,unsigned int>(aij,ptr));
        }

        std::list<unsigned int> elms_at = NONEMPTY[i-1];
        for (std::list<unsigned int>::iterator it=elms_at.begin(); it!=elms_at.end(); ++it) {
            unsigned int ajj = *it;
            std::pair<boost::unordered_multimap<unsigned int, unsigned int>::iterator,boost::unordered_multimap<unsigned int, unsigned int>::iterator> ret = Q.equal_range(ajj);

            for (boost::unordered_multimap<unsigned int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
                const unsigned int &ptr = it->second;
                QUEUE.push_back(ptr);
            }
        }
        Q.clear();
    }

    for (std::deque<unsigned int>::iterator itr=QUEUE.begin(); itr!=QUEUE.end(); ++itr) {
        unsigned int ptr = *itr;
        s_ptrs.push_back(ptr);
    }
}


int isomorphism(const aw::Tree &treei, const aw::Tree &treek, const aw::TreetaxaMap &g_nmapi, const aw::TreetaxaMap &g_nmapk) {

    std::vector<aw::Tree> g_trees;
    g_trees.push_back(treei);
    g_trees.push_back(treek);

    std::vector<aw::TreetaxaMap> g_nmaps;
    g_nmaps.push_back(g_nmapi);
    g_nmaps.push_back(g_nmapk);

    if(g_nmaps[0].tot_leaves() != g_nmaps[1].tot_leaves()) {        
        return 1;             //Not Isomorphic when total number of leaves are different
    }
    
    if(g_nmaps[0].unq_leaves() != g_nmaps[1].unq_leaves()) {
        return 1;             //Not Isomorphic when unique leaves are different
    }
    
    //check if each global id has smae number of duplicates
    std::set<unsigned int> gids;
    g_nmaps[0].unq_gids(gids);
    for (std::set<unsigned int>::iterator it=gids.begin(); it!=gids.end(); ++it) {
        unsigned int ggid =  *it;
        int copies0 = g_nmaps[0].ids_count(ggid);
        int copies1 = g_nmaps[1].ids_count(ggid);
        if(copies0 != copies1)
            return 1;
    }

    //computing heights for trees..
    for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
        int hgt = 0;
        int maxhgt = 0;

        for (aw::Tree::iterator_dfs m=g_trees[i].begin_dfs(),mEE=g_trees[i].end_dfs(); m!=mEE; ++m) {
            if(m.idx == g_trees[i].root) {
                hgt = 0;
                continue;
            }
            switch (m.direction) {
                case aw::PREORDER: {
                    hgt++;
                    if(g_trees[i].is_leaf(m.idx))
                        if(maxhgt < hgt)
                            maxhgt = hgt;
                } break;
                case aw::POSTORDER: {
                    if(g_trees[i].is_leaf(m.idx))
                        continue;
                    hgt--;
                } break;
                case aw::INORDER: {
                    if(g_trees[i].is_leaf(m.idx))
                        continue;
                    hgt--;
                } break;
                default: break;
            }
        }
        g_trees[i].height = maxhgt;
    }

    if(g_trees[0].height != g_trees[1].height) {
        return 1;        //Not Isomorphic when heights are different
    }
    
    //storing nodes for levels...
    for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
        int hgt = 0;
        int level;
        for (aw::Tree::iterator_dfs m=g_trees[i].begin_dfs(),mEE=g_trees[i].end_dfs(); m!=mEE; ++m) {
            if(m.idx == g_trees[i].root) {
                hgt = 0;
                level = g_trees[i].height;
                g_trees[i].level2id.insert(std::pair<int,unsigned int>(level,m.idx));
                g_trees[i].set_level(m.idx,level);
                continue;
            }
            switch (m.direction) {
                case aw::PREORDER: {
                    hgt++;
                    level = g_trees[i].height - hgt;
                    g_trees[i].level2id.insert(std::pair<int,unsigned int>(level,m.idx));
                    g_trees[i].set_level(m.idx,level);
                } break;
                case aw::POSTORDER: {
                    if(g_trees[i].is_leaf(m.idx))
                        continue;
                    hgt--;
                } break;
                case aw::INORDER: {
                    if(g_trees[i].is_leaf(m.idx))
                        continue;
                    hgt--;
                } break;
                default: break;
            }
        }
    }

    for(int h=0,hEE=g_trees[0].height; h<=hEE; ++h) {
        std::vector<std::vector<unsigned int> > listsB;
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            std::vector<unsigned int> ptr;
            std::vector<unsigned int> lists;
            if(h==0) {  //indeces will bboost::unordered_multimap<int, unsigned int>::iteratore leaf labels (integer lebels) at level 0
                std::pair<boost::unordered_multimap<int, unsigned int>::iterator,boost::unordered_multimap<int, unsigned int>::iterator> ret = g_trees[i].level2id.equal_range(h);
                for (boost::unordered_multimap<int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
                    const unsigned int &id = it->second;
                    g_trees[i].update_index(id,g_nmaps[i].gid(id));
                    ptr.push_back(lists.size());
                    lists.push_back(1);
                    lists.push_back(g_nmaps[i].gid(id));
                }
            }
            else {
                std::pair<boost::unordered_multimap<int, unsigned int>::iterator,boost::unordered_multimap<int, unsigned int>::iterator> ret = g_trees[i].level2id.equal_range(h);
                for (boost::unordered_multimap<int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
                    const unsigned int &id = it->second;
                    std::list<unsigned int> indices;
                    if(g_trees[i].is_leaf(id))
                        indices.push_back(g_nmaps[i].gid(id));
                    else {
                        BOOST_FOREACH(const unsigned int &u, g_trees[i].adjacent(id)) {
                            if (g_trees[i].return_level(u) == h-1)
                                indices.push_back(g_trees[i].return_index(u));
                        }
                        indices.sort();
                    }
                    int len_ind = indices.size();
                    ptr.push_back(lists.size());
                    lists.push_back(len_ind);
                    for (std::list<unsigned int>::iterator it = indices.begin(); it != indices.end(); it++)
                        lists.push_back(*it);
                }
            }

            //sort the tuples at this level
            std::vector<unsigned int> sort_ptrs;
            sort_laxi(ptr,lists,sort_ptrs);


            std::vector<unsigned int> new_lists;
            BOOST_FOREACH(const unsigned int &c, sort_ptrs) {
                int tot_elms = lists[c];
                for(int j = c; j<=c+tot_elms; ++j)
                    new_lists.push_back(lists[j]);
            }
            listsB.push_back(new_lists);

            //update indices of tree-nodes if h != 0
            if(h != 0) {
                boost::unordered_map<unsigned int,unsigned int> sptr2rank; //sorted_pointer to ranks

                std::vector<unsigned int> initt;
                unsigned int f = sort_ptrs[0];
                int elms = lists[f];
                for(int fk = f+1; fk<=f+elms; ++fk) {
                    initt.push_back(lists[fk]);
                }

                unsigned rnk = 0;
                sptr2rank.insert(std::pair<unsigned int,unsigned int>(f,rnk)); //first has first rank
                for(int fk=1; fk<sort_ptrs.size(); ++fk) {
                    f = sort_ptrs[fk];
                    elms = lists[f];

                    if(elms!=initt.size())  {
                        sptr2rank.insert(std::pair<unsigned int,unsigned int>(f,++rnk));
                    }
                    else {
                        int flag = 0;
                        int m = 0;
                        for(int fk = f+1; fk<=f+elms; ++fk) {
                            if(lists[fk]!=initt[m++]) {
                                flag = 1;
                                break;
                        }    }
                        if(flag==0) {
                            sptr2rank.insert(std::pair<unsigned int,unsigned int>(f,rnk));
                            continue;
                        }
                        else {
                            sptr2rank.insert(std::pair<unsigned int,unsigned int>(f,++rnk));
                        }
                    }

                    //update initt
                    initt.clear();
                    for(int fk = f+1; fk<=f+elms; ++fk)
                        initt.push_back(lists[fk]);
                }

                //update indices now
                std::pair<boost::unordered_multimap<int, unsigned int>::iterator,boost::unordered_multimap<int, unsigned int>::iterator> ret = g_trees[i].level2id.equal_range(h);
                int pp = 0;
                //std::cout<<"\n Tree:"<<i;
                for (boost::unordered_multimap<int, unsigned int>::iterator it=ret.first; it!=ret.second; ++it) {
                    const unsigned int &id = it->second;
                    unsigned int rnk = sptr2rank[ptr[pp]];
                    g_trees[i].update_index(id,rnk);
                    //std::cout<<"<"<<id<<" "<<rnk<<">";
                    pp++;
                }


            }
        }

        if (h==g_trees[0].height && g_trees[0].return_index(0)!=g_trees[1].return_index(0)) {
            std::cout<<g_trees[0].return_index(0);
            //std::cout<<" "<<g_trees[1].return_index(0);
            //MSG("\nNot Isomorphic !!!"<<h);
            return 1;
        }

        //compare sorted list of tuples
        int ln = listsB[1].size();
        std::vector<unsigned int> ls0 = listsB[0];
        std::vector<unsigned int> ls1 = listsB[1];
        for(int i = 0; i<ln; ++i) {
            if(ls0[i]!= ls1[i]) {
                //MSG("\nNot Isomorphic !!!"<<h);
                return 1;
            }
        }
    }

    return 0;

}

#endif /*ISO_C_*/


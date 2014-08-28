//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		CHashFunc: Universal hash functions
//      ***** Updated 2013 : Ruchi Chaudhary (UFL) for Mul-tree compression ****
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

#ifndef HASHFUNC_HH
#define HASHFUNC_HH

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <boost/random.hpp>

#include "global.h"

typedef struct {
	unsigned long long hv1;
	unsigned long long hv2;
} HV_STRUCT_T;	

class CHashFunc {
	
	unsigned long long    _m1;   // prime number1 for hash function1
	unsigned long long    _m2;   // prime number1 for hash function2
	unsigned int          _t;    // number of trees
	unsigned int          _n;    // number of taxa
	unsigned long long *  _a1;   // random numbers for hash function1
	unsigned long long *  _a2;   // random numbers for hash function2
	unsigned int	      _c;    // double collision factor: constant for c*t*n of hash function2;
	
public:
	CHashFunc(unsigned int t, unsigned int n, unsigned int c) : _m1(0), _m2(0), _t(t), _n(n), _a1(NULL), _a2(NULL), _c(c)
        {
            UHashfunc_init(t, n, c);
        }
        
	int UHashfunc_init(unsigned int t, unsigned int n, unsigned int c) {
            // Init member variables
            _t = t;  _n = n;  _c = c;
            
            // Get the prime number which is larger than t*n
            unsigned long long top = _t*_n;
            unsigned long long p = 0;
            unsigned int mul = 1;
            do {
                unsigned from = 100 * mul;
                p = GetPrime(top, from);
                ++mul;
            } while (p == 0);
            _m1 = p;   	           

            unsigned long long top2 = _c*_t*_n;
            unsigned long long p2 = 0;
            mul = 1;
            do {
                unsigned from = 100 * mul;
                p2 = GetPrime(top2, from);
                ++mul;
            } while (p2 == 0);

            _m2 = p2;           

            _a1 = new unsigned long long[_n];
            _a2 = new unsigned long long[_n];

            for (unsigned int i=0; i<_n; ++i) {
                boost::uniform_int<> range1(0,_m1);
                boost::variate_generator<boost::mt19937&, boost::uniform_int<> > getRnd(rng, range1);
                _a1[i] = getRnd();                

                boost::uniform_int<> range2(0,_m2);
                boost::variate_generator<boost::mt19937&, boost::uniform_int<> > getRnd2(rng, range2);
                _a2[i] = getRnd2();                
            }          
        }
        
        ~CHashFunc() {
            delete[] _a1;  delete[] _a2;
        }	
	
	unsigned long long GetPrime(unsigned long long topNum, unsigned from){
            unsigned long long primeNum=0;
            unsigned long long candidate=0;

            if (topNum <= 100)  candidate = 2;
            else   candidate = topNum;

            while (candidate <= topNum+from) {
                unsigned long long trialDivisor = 2;
                int prime = 1;

                while (trialDivisor * trialDivisor <= candidate) {
                        if (candidate % trialDivisor == 0) {
                                prime = 0;
                                break;
                        }
                        trialDivisor++;
                }
                if (prime) primeNum = candidate;
                candidate++;
            }
            return primeNum;
        }
	
	//Implicit bp
	unsigned long long getA1(unsigned idx) { return (_a1[idx]); }
	unsigned long long getA2(unsigned idx) { return (_a2[idx]); }
	unsigned long long getM1() { return _m1; }
	unsigned long long getM2() { return _m2; }

        CHashFunc() {}
};

#endif

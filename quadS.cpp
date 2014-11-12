//
//  qs.cpp
//  avalgproj1
//
//  Created by Tobias Wikström on 12/11/14.
//  Copyright (c) 2014 Tobias Wikström. All rights reserved.
//

#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <mpfr.h>
//#include "gmpxx.h"
using namespace std;


class quadS {
    
private:
    
    vector<int> basePrimes;
    mpz_t N;
    mpz_t rootN;
    mpz_t rootN2;
    mpz_t qValLastX;
    mpz_t lastValueForXOffset;
    
public:
    
    
    quadS() {}
    
    ~quadS() {}
    
    quadS(int s_N) {
        
    }
    
    
    quadS(char * s_N) {
        mpz_t N;
        mpz_set_str(N, s_N, 10);
        init(N);
    }
    
    void init(mpz_t n) {
        mpz_init(N);
        mpz_init(rootN);
        mpz_init(rootN2);
        mpz_init(lastValueForXOffset);
        mpz_set(N, n);
        
        mpz_sqrt(rootN, n);
        
        mpz_add(rootN2, rootN, rootN);
        mpz_t temp;
        mpz_init(temp);
        mpz_mul(temp, rootN, rootN);
        mpz_sub(qValLastX, temp, n);
        
        
    }
    
    /* 1. Generate factor base */
    
    void generateFactorBase(mpz_t B, mpz_t N) {
        mpfr_t fN, lnN, lnlnN;
        mpfr_init(fN), mpfr_init(lnN), mpfr_init(lnlnN);
        
        mpfr_set_z(fN, N, MPFR_RNDU);
        mpfr_log(lnN, fN, MPFR_RNDU);
        mpfr_log(lnlnN, lnN, MPFR_RNDU);
        
        mpfr_mul(fN, lnN, lnlnN, MPFR_RNDU);
        mpfr_sqrt(fN, fN, MPFR_RNDU);
        mpfr_div_ui(fN, fN, 2, MPFR_RNDU);
        mpfr_exp(fN, fN, MPFR_RNDU);
        
        mpfr_get_z(B, fN, MPFR_RNDU);
        
        mpfr_clears(fN, lnN, lnlnN, NULL);
        
    }
    
    void printMembers() {
        char * s;
        s = mpz_get_str(NULL, 10, rootN);
        cout << "Sqrt of N:  " <<  s << endl;
        s = mpz_get_str(NULL, 10, rootN2);
        cout << "2 * Sqrt of N:  " <<  s << endl;
        s = mpz_get_str(NULL, 10, qValLastX);
        cout << "qVal last X:  " <<  s << endl;
        s = mpz_get_str(NULL, 10, lastValueForXOffset);
        cout << "lastValueForXOffset:  " <<  s << endl;
        
    }
    
    void print(mpz_t a) {
        char * s;
        s = mpz_get_str(NULL, 10, a);
        cout <<  s << endl;

    }
    
    void gcd(mpz_t & ret, mpz_t a, mpz_t b) {
		mpz_t c;
        mpz_init(c);
        mpz_t zero;
        mpz_init(zero);
		while (mpz_cmp(a, zero) != 0) {
            mpz_set(c, a);
            mpz_mod(a, b, a);
            mpz_set(b, c);
		}
		mpz_set(ret, b);
	}
    
};

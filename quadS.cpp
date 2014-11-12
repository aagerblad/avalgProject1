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
using namespace std;

#define GET_BIT_AT(index) ((numbers[index>>3] & (1<<(index&7))) >> (index&7))
#define SET_BIT_AT(index) (numbers[index>>3] |= (1 << (index&7)))
#define CLEAR_BIT_AT(index) (numbers[index>>3] &= ~(1 << (index&7)))


class quadS {
    
    
private:
    
    typedef struct modular_root {
        unsigned long root1;
        unsigned long root2;
    } modular_root_t;
    
    // Quadtratic sieve members
    vector<int64_t> basePrimes;
    vector<modular_root> modular_roots;
    mpz_t N;
    mpz_t rootN;
    mpz_t rootN2;
    mpz_t qValLastX;
    mpz_t lastValueForXOffset;
    int64_t quadraticPrimesFound;
    
    // Eratosthenes members
    char *numbers;
    int64_t base_ref;
    
public:
    
    
    quadS() {}
    
    ~quadS() {
//        delete[] modular_roots;
    }
    
    quadS(mpz_t N) {
        init(N);
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
        mpz_init(qValLastX);
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
    
    void generateModularRoots() {
//        modular_roots = new modular_root[quadraticPrimesFound];
        mpz_t tmp, r1, r2;
        mpz_init(tmp);
        mpz_init(r1);
        mpz_init(r2);
        modular_roots.resize(quadraticPrimesFound);
        
        for (int i = 0; i < quadraticPrimesFound; i++) {
            mpz_set_ui(tmp, (unsigned long) basePrimes[i]);
            mpz_sqrtrem(r1, N, tmp); /* calculate the modular root */
            mpz_neg(r2, r1); /* -q mod n */
            mpz_mod(r2, r2, tmp);

            modular_roots[i].root1 = mpz_get_ui(r1);
            modular_roots[i].root2 = mpz_get_ui(r2);
        }
        cout << modular_roots.size() << endl;
        mpz_clear(tmp);
        mpz_clear(r1);
        mpz_clear(r2);
    }
    

    /************************************************/
    /*              Eratosthenes sieve              */
    /************************************************/
    
    
    //Can handle up to base 2^31 * 8 * 8
    int64_t sieve_primes_up_to(int64_t base) {
        int64_t num_primes = 0;
        int64_t i;
        
        //Find primes, base included
        base++;
        
        numbers = (char*) calloc(base / 64 + 1, sizeof(uint64_t));
        base_ref = base;
        
        for (i = 0; i < base; i++)
            SET_BIT_AT(i);
        
        int64_t p = 2;
        num_primes++;
        int64_t offset;
        
        while (1) {
            offset = 2 * p;
            
            while (offset < base) {
                CLEAR_BIT_AT(offset);
                offset += p;
            }
            
            offset = p + 1;
            while (offset < base && (GET_BIT_AT(offset) == 0)) {
                offset++;
            }
            
            if (offset == base)
                break;
            
            p = offset;
            num_primes++;
        }
        
        return num_primes;
    }
    
    /* Array must be malloc'ed already
     * fills the array with the first @num_primes primes calculated with sieve_primes
     */
    void fill_primes(int64_t *primes_array) {
        int64_t j, i;
        
        for (j = 0, i = 2; i < base_ref; i++) {
            if (GET_BIT_AT(i) == 1) {
                primes_array[j] = i;
                j++;
            }
        }
        
        free(numbers);
    }
    
    /* Fill the array with only primes where n is a quadratic residue: x² = n (mod p) */
    int64_t fill_primes_with_quadratic_residue(mpz_t n) {
        int64_t j, i;
        mpz_t b;
        mpz_init(b);
        
        basePrimes.push_back(2);
        for (j = 1, i = 3; i < base_ref; i++) {
            mpz_set_ui(b, (unsigned long) i);
            if ((GET_BIT_AT(i)) == 1 && mpz_jacobi(n, b) == 1) {
//                basePrimes[j] = i;
                basePrimes.push_back(i);
                j++;
            }
        }
        
        free(numbers);
        quadraticPrimesFound = j;
        return j;
    }
    
    /************************************************/
    /*               Helper functions               */
    /************************************************/
    
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
    
    int64_t toInt(mpz_t v) {
        return mpz_get_ui(v);
    }
    
    vector<int64_t> getPrimes() {
        return basePrimes;
    }
    
    vector<modular_root> getModularRoots() {
        return modular_roots;
    }
    
    void printModularRoots() {
        for (modular_root r : modular_roots) {
            cout << r.root1 << " " << r.root2 << endl;;
        }
    }
    
};

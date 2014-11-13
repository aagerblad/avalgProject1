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
    
    typedef struct {
        mpz_t *MATRIX;
        mpz_t *IDENTITY;
        uint64_t rows;
        uint64_t cols;
        uint64_t next_free_row;
    } matrix_t;
    
    typedef struct {
        mpz_t value_x;
        mpz_t value_x_squared;
        
        mpz_t factors_vect;
        
    } smooth_number_t;
    
    // Quadtratic sieve members
    vector<int64_t> basePrimes;
    vector<modular_root> modular_roots;
    mpz_t N;
    mpz_t rootN;
    mpz_t rootN2;
    mpz_t qValLastX;
    mpz_t lastValueForXOffset;
    int64_t quadraticPrimesFound;
    uint64_t smoothNmbrsFound = 0;
    
    matrix_t matrix;
    int vec_offset;
    
    smooth_number_t *smooth_numbers;
    
    // Eratosthenes members
    char *numbers;
    int64_t base_ref;
    
public:
    
    
    quadS() {}
    
    ~quadS() {
//        delete[] modular_roots;
        delete[] matrix.MATRIX;
        delete[] matrix.IDENTITY;
        delete[] numbers;
        delete[] smooth_numbers;
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
    
    /************************************************/
    /*               Matrix functions               */
    /************************************************/
    
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
    
#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
    
    /* allocates space for m rows * n columns matrix for MATRIX and IDENTITY */
    void matrix_init(uint64_t m, uint64_t n) {
        m = m + vec_offset;
        matrix.MATRIX = new mpz_t[m]; //(m, sizeof(mpz_t));
        matrix.IDENTITY = new mpz_t[m]; //(mpz_t *) calloc(m, sizeof(mpz_t));
        matrix.rows = m;
        matrix.cols = n;
        matrix.next_free_row = 0;
        mpz_init2(tmp_matrix_row, quadraticPrimesFound);
    }
    
    void matrix_push_row(mpz_t row) {
        mpz_init(matrix.MATRIX[matrix.next_free_row]);
        mpz_set(matrix.MATRIX[matrix.next_free_row], row);
        
        mpz_init2(matrix.IDENTITY[matrix.next_free_row], matrix.cols); /* initializes a n bit vector all set to 0 */
        mpz_setbit(matrix.IDENTITY[matrix.next_free_row], matrix.next_free_row); /* set the next_free_row bit to 1 */
        matrix.next_free_row++;
    }
    
    void matrix_print_matrix() {
        uint64_t i, j;
        
        printf("\nMATRIX\n");
        for (i = 0; i < matrix.rows; i++) {
            printf("[");
            for (j = 0; j < matrix.cols; j++) {
                printf("%2d", mpz_tstbit(matrix.MATRIX[i], j));
            }
            printf(" ]\n");
        }
        printf("\n");
    }
    
    void matrix_print_identity() {
        uint64_t i, j;
        
        printf("\nIDENTITY\n");
        for (i = 0; i < matrix.rows; i++) {
            printf("[");
            for (j = 0; j < matrix.cols; j++) {
                printf("%2d", mpz_tstbit(matrix.IDENTITY[i], j));
            }
            printf(" ]\n");
        }
        printf("\n");
    }
    
    void matrix_free() {
        free(matrix.MATRIX);
        free(matrix.IDENTITY);
    }
    
    /* performs a Gauss elimination on matrix->MATRIX, result (linear dependence) will be in the matrix->IDENTITY */
    void gauss_elimination() {
        printf("\nPerforming Gauss elimination..\n");
        mpz_t *m = matrix.MATRIX;
        mpz_t *I = matrix.IDENTITY;
        
        uint64_t col, row, next_row, next_pivot;
        for (next_row = 0, col = 0; col < min(matrix.cols, matrix.rows); col++) /* for all rows*/
        {
            next_pivot = -1;
            for (row = next_row; row < matrix.rows; row++) /* search for the next pivot*/
            {
                if (mpz_tstbit(m[row], col)) {
                    next_pivot = row; /* row contains the next pivot */
                    next_row++;
                    break;
                }
            }
            
            if (next_pivot == -1)
                continue;
            
            if (next_pivot != next_row - 1) /* current row is not the pivot, switch rows */
            {
                mpz_swap(m[next_pivot], m[next_row - 1]);
                mpz_swap(I[next_pivot], I[next_row - 1]);
            }
            
            for (row = next_row; row < matrix.rows; row++) {
                if (mpz_tstbit(m[row], col)) {
                    mpz_xor(m[row], m[row], m[next_row - 1]); /* XOR the rows to eliminate the 1 in position (row, next_row-1)*/
                    mpz_xor(I[row], I[row], I[next_row - 1]);
                }
            }
        }
    }
    
    /* does not check for bounds, the caller must */
    void get_matrix_row(mpz_t rop, uint64_t row_index) {
        mpz_set(rop, matrix.MATRIX[row_index]);
    }
    
    void get_identity_row(mpz_t rop, uint64_t row_index) {
        mpz_set(rop, matrix.IDENTITY[row_index]);
    }
    
    
    /************************************************/
    /*          Quadratic sieve functions           */
    /************************************************/
    
    
    
    /* 1. Generate factor base */
    // TODO kanske ha en tabell istället? Äh, tar inte så långt tid tror jag.
    
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
    
    /* 2. Generate modular roots (solve equation) */
    
//    void CRT(mpz_t res, mpz_t N, mpz_t p){
//        
//        mpz_t M;
//        mpz_set_ui(M, 1);
//        mpz_mul(M, M, p);
//        vector< long long > m, s;
//        for(int i=0; i<int(mods.size()); i++){
//            m.push_back(M/mods[i]);
//            long long temp=m[i]%mods[i];
//            long long k=0;
//            /* if there is a possibility of k being very big, then prime factorize m[i],
//             * find modular inverse of 'temp' of each of the factors
//             * 'k' equals to the multiplication ( modular mods[i] ) of modular inverses
//             */
//            while(true){
//                if((k*temp)%mods[i]==1) break;
//                k++;
//            }
//            s.push_back(k);
//        }
//        long long ret=0;
//        for(int i=0; i<int(s.size()); i++) {
//            ret+=( (m[i]*s[i])%M *r[i] )%M;
//            if(ret>=M) ret-=M;
//        }
//    }
    
    void legendre(mpz_t r1, mpz_t N, mpz_t p) {
        mpz_t a;
        mpz_init(a);
        mpz_set_ui(r1, 1);
        
        mpz_mod(a, N, p);
        unsigned long pu;
        pu = mpz_get_ui(p);
        long power = (pu-1)/2;
        
        while (power > 0) {
            if (power % 2 == 1) {
//                unsigned int res = 1;
                mpz_mul(r1, r1, a);
                mpz_mod(r1, r1, p);
//                mpz_set_ui
            }
            
            mpz_mul(a, a, a);
            mpz_mod(a, a, p);
            power = power / 2;
        }
        
        mpz_t tmp;
        mpz_init(tmp);
        mpz_sub(tmp, r1, p);
        if (mpz_cmp_ui(tmp, -1) == 0) {
            mpz_sub(r1, r1, p);
        }
        
        
    }
    
    int shanksAndTonelli(mpz_t q, const mpz_t n, const mpz_t p) {
        mpz_t w, n_inv, y;
        unsigned int i, s;

        if(mpz_divisible_p(n, p)) {        
            mpz_set_ui(q, 0);
            return 1;
        }
   
        if(mpz_tstbit(p, 1) == 1) {
            mpz_set(q, p);
            mpz_add_ui(q, q, 1);
            mpz_fdiv_q_2exp(q, q, 2);
            mpz_powm(q, n, q, p);
            return 1;
        }
        mpz_init(y);
        mpz_init(w);
        mpz_init(n_inv);
        
        mpz_set(q, p);
        mpz_sub_ui(q, q, 1);
        s = 0;
        while(mpz_tstbit(q, s) == 0) s++;
        mpz_fdiv_q_2exp(q, q, s);
        mpz_set_ui(w, 2);
        while(mpz_legendre(w, p) != -1)
            mpz_add_ui(w, w, 1);
        mpz_powm(w, w, q, p);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 1);
        mpz_powm(q, n, q, p);
        mpz_invert(n_inv, n, p);
        for(;;) {
            mpz_powm_ui(y, q, 2, p);
            mpz_mul(y, y, n_inv);
            mpz_mod(y, y, p);
            i = 0;
            while(mpz_cmp_ui(y, 1) != 0) {
                i++;
                mpz_powm_ui(y, y, 2, p);
            }
            if(i == 0) {
                return 1;
            }
            if(s-i == 1) {
                mpz_mul(q, q, w);
            } else {
                mpz_powm_ui(y, w, 1 << (s-i-1), p);
                mpz_mul(q, q, y);
            }
            mpz_mod(q, q, p);
        }
        
        mpz_clear(w); mpz_clear(n_inv); mpz_clear(y);
        return 0;
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
            
            shanksAndTonelli(r1, N, tmp);
            mpz_neg(r2, r1); /* -q mod n */
            mpz_mod(r2, r2, tmp);
            
            cout << "Base P: " << basePrimes[i];
            
            modular_roots[i].root1 = mpz_get_ui(r1);
            modular_roots[i].root2 = mpz_get_ui(r2);
            
            cout << " <--> " << modular_roots[i].root1 << " and " << modular_roots[i].root2 << endl;
        }
        cout << "Modular roots: " << modular_roots.size() << endl;
        mpz_clear(tmp);
        mpz_clear(r1);
        mpz_clear(r2);
    }
    
    /* 3. Sieving */
    
    // returns the index of the first element to start sieving from
    // (first multiple of root that is directly greater than start)
    // res = p*t + root >= start */
    void get_sieving_start_index(mpz_t res, mpz_t start, mpz_t p,
                                 unsigned long root) {
        mpz_t q, r;
        mpz_init(q);
        mpz_init(r);
        
        mpz_sub_ui(start, start, root);
        mpz_fdiv_qr(q, r, start, p);
        
        if (mpz_cmp_ui(r, 0) != 0)
            mpz_add_ui(q, q, 1);
        
        mpz_mul(q, q, p); /* next element p*q+root that is directly >= start */
        mpz_add_ui(q, q, root);
        mpz_set(res, q);
        mpz_clear(q);
        mpz_clear(r);
    }
    
    mpz_t tmp_matrix_row;
    void save_smooth_number(smooth_number_t n) {
        mpz_clear(tmp_matrix_row); /* tmp_matrix_row must be initialized already */
        mpz_init2(tmp_matrix_row, quadraticPrimesFound); /* init a vector of *exactly* nb_qr_primes bits */
        
        if (smoothNmbrsFound > quadraticPrimesFound + vec_offset - 1) /* if we have sufficient smooth numbers, skip saving */
            return;
        
        smooth_number_t tmp;
        mpz_init(tmp.value_x);
        mpz_init(tmp.value_x_squared);
        
        mpz_set(tmp.value_x, n.value_x);
        mpz_pow_ui(tmp.value_x_squared, n.value_x, 2);
        mpz_sub(tmp.value_x_squared, tmp.value_x_squared, N); /* x²-N. Saving this will enable us to not go through exponents
                                                               * and reconstruct the original number */
        
        smooth_numbers[smoothNmbrsFound++] = tmp;

        
        /* the coefficient vector in GF2 has already been constructed */
        matrix_push_row(n.factors_vect);
    }
    
    void sieve() {
        mpz_t x, sieving_index, next_sieving_index;
        unsigned long SIEVING_STEP = 50000; /* we sieve for 50000 elements at each loop */
        uint64_t p_pow;
        smooth_number_t *x_squared;
        
        x_squared = new smooth_number_t[SIEVING_STEP];
        
        smooth_numbers = new smooth_number_t[quadraticPrimesFound];
        
        mpz_init_set(x, rootN);
        mpz_init_set(sieving_index, x);
        mpz_init_set(next_sieving_index, x);
        
        // Init
        
        mpz_t p;
        mpz_init(p);
        mpz_t str;
        mpz_init_set(str, sieving_index);
        
        for (int i = 0; i < SIEVING_STEP; i++) {
            mpz_init(x_squared[i].value_x);
            mpz_init(x_squared[i].value_x_squared);
            
            mpz_init2(x_squared[i].factors_vect, quadraticPrimesFound);
            mpz_add_ui(x, x, 1);
        }
        
        int nb_smooth_per_round = 0;
//        char s[512];
        
        cout << "sieving..." << endl;
        
        while (smoothNmbrsFound < quadraticPrimesFound + vec_offset) {
            nb_smooth_per_round = 0;
            mpz_set(x, next_sieving_index); /* sieve numbers from sieving_index to sieving_index + sieving_step */
            mpz_set(sieving_index, next_sieving_index);
            
//            fflush(stdout);
//            cout << "hej" << endl;
            for (int i = 0; i < SIEVING_STEP; i++) {
                mpz_set(x_squared[i].value_x, x);
                
                mpz_pow_ui(x_squared[i].value_x_squared, x, 2); /* calculate value_x_squared <- x²-n */
                mpz_sub(x_squared[i].value_x_squared, x_squared[i].value_x_squared,
                        N);
                
                mpz_clear(x_squared[i].factors_vect);
                mpz_init2(x_squared[i].factors_vect, quadraticPrimesFound); /* reconstruct a new fresh 0ed vector of size nb_qr_primes bits */
                
                mpz_add_ui(x, x, 1);
            }
            mpz_set(next_sieving_index, x);
            
            for (int i = 0; i < quadraticPrimesFound; i++) {
                mpz_set_ui(p, (unsigned long) basePrimes[i]);
                mpz_set(x, sieving_index);
                
                /* get the first multiple of p that is directly larger that sieving_index
                 * Quadratic SIEVING: all elements from this number and in positions multiples of root1 and root2
                 * are also multiples of p */
                get_sieving_start_index(x, x, p, modular_roots[i].root1);
                mpz_set(str, x);
                mpz_sub(x, x, sieving_index); /* x contains index of first number that is divisible by p */
                
                for (unsigned long j = mpz_get_ui(x); j < SIEVING_STEP; j += basePrimes[i]) {
                    p_pow = mpz_remove(x_squared[j].value_x_squared,
                                       x_squared[j].value_x_squared, p); /* eliminate all factors of p */
                    
                    if (p_pow & 1) /* mark bit if odd power of p exists in this x_squared[j] */
                    {
                        mpz_setbit(x_squared[j].factors_vect, i);
                    }
//                    print(x_squared[j].value_x); // TODO remove
                    
                    if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0) {
                        cout << "saving smooth" << endl;
                        save_smooth_number(x_squared[j]);
                        nb_smooth_per_round++;
                    }
                    /* sieve next element located p steps from here */
                }
                
                /* same goes for root2 */
                if (modular_roots[i].root2 == modular_roots[i].root1)
                    continue;
                
                mpz_set(x, sieving_index);
                
                get_sieving_start_index(x, x, p, modular_roots[i].root2);
                mpz_set(str, x);
                mpz_sub(x, x, sieving_index);
                
                for (unsigned long j = mpz_get_ui(x); j < SIEVING_STEP; j += basePrimes[i]) {
                    p_pow = mpz_remove(x_squared[j].value_x_squared,
                                       x_squared[j].value_x_squared, p);
                    
                    if (p_pow & 1) {
                        mpz_setbit(x_squared[j].factors_vect, i);
                    }
                    
                    if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0) {
                        save_smooth_number(x_squared[j]);
                        nb_smooth_per_round++;
                    }
                }
            }
        }
    }
    
    /* Factor */
    
    void factor() {
        uint64_t row_index = quadraticPrimesFound + vec_offset - 1; /* last row in the matrix */
        int nb_linear_relations = 0;
        mpz_t linear_relation_z, solution_z;
        mpz_init(linear_relation_z);
        mpz_init(solution_z);
        
        get_matrix_row(linear_relation_z, row_index--); /* get the last few rows in the Gauss eliminated matrix*/
        while (mpz_cmp_ui(linear_relation_z, 0) == 0) {
            nb_linear_relations++;
            get_matrix_row(linear_relation_z, row_index--);
        }
        
        printf("\nFactorizing..\n");
        mpz_t solution_X, solution_Y;
        mpz_init(solution_X);
        mpz_init(solution_Y);
        
        /* we start testing from the first linear relation encountered in the matrix */
        for (int j = nb_linear_relations; j > 0; j--) {
            printf("Trying %d..\n", nb_linear_relations - j + 1);
            mpz_set_ui(solution_X, 1);
            mpz_set_ui(solution_Y, 1);
            
            get_identity_row(solution_z,
                             quadraticPrimesFound + vec_offset - j + 1);
            
            for (int i = 0; i < quadraticPrimesFound; i++) {
                if (mpz_tstbit(solution_z, i)) {
                    mpz_mul(solution_X, solution_X, smooth_numbers[i].value_x);
                    mpz_mod(solution_X, solution_X, N); /* reduce x to modulo N */
                    
                    mpz_mul(solution_Y, solution_Y,
                            smooth_numbers[i].value_x_squared);
                    /*TODO: handling huge stuff here, there is no modulo N like in the solution_X case!
                     * eliminate squares as long as you go*/
                }
            }
            
            mpz_sqrt(solution_Y, solution_Y);
            mpz_mod(solution_Y, solution_Y, N); /* y = sqrt(MUL(xi²-n)) mod N */
            
            mpz_sub(solution_X, solution_X, solution_Y);
            
            mpz_gcd(solution_X, solution_X, N);
            
            if (mpz_cmp(solution_X, N) != 0 && mpz_cmp_ui(solution_X, 1) != 0) /* factor can be 1 or N, try another relation */
                break;
        }
        
        unsigned long hej = mpz_get_ui(solution_X);
        cout << "Warning: Division by " << hej << endl;
        mpz_cdiv_q(solution_Y, N, solution_X);
        
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
        
        numbers = new char[base / 64 + 1];
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

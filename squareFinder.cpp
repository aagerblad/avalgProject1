#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "primenumbers.cpp"

using namespace std;

class squareFinder
{
public:
	squareFinder();
	~squareFinder();

	/* data */
	vector<mpz_class> getSquares(mpz_class N) {



		bitset<B> matrix[B];

		cout << matrix[1] << endl;


	}	


	mpz_class generateBaseLimit(mpz_class N) {

        mpfr_t fN, lnN, lnlnN;
        mpfr_init(fN), mpfr_init(lnN), mpfr_init(lnlnN);

        mpfr_set_z(fN, N, MPFR_RNDU);
        mpfr_log(lnN, fN, MPFR_RNDU);
        mpfr_log(lnlnN, lnN, MPFR_RNDU);
        
        mpfr_mul(fN, lnN, lnlnN, MPFR_RNDU);
        mpfr_sqrt(fN, fN, MPFR_RNDU);
        mpfr_div_ui(fN, fN, 2, MPFR_RNDU);
        mpfr_exp(fN, fN, MPFR_RNDU);
        
        mpz_t B; mpz_t.init(B);
        mpfr_get_z(B, fN, MPFR_RNDU);

        mpfr_clears(fN, lnN, lnlnN, NULL);

        return mpz_class(B);
        
    }
};
#include <iostream>
#include <queue>
#include <gmpxx.h>
#include "quadS.cpp"
#include "miller_rabin.h"

#include "primenumbers.cpp" // warning: This include will allocate a bunch of memory
using namespace std;

void printis(mpz_t a) {
    char * s;
    s = mpz_get_str(NULL, 10, a);
    cout <<  s;
    
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


mpz_t xSquared, a, n;

//rop <- (rop² + a) mod n
void f(mpz_t rop)
{
	mpz_pow_ui(xSquared, rop, 2);
	mpz_add(xSquared, xSquared, a);
	mpz_mod(rop, xSquared, n);
}

void runPollard(mpz_t nisse, mpz_t x, mpz_t y)
{
    
	mpz_t g, gMul, s, u, v, uMinusV;
	mpz_t nMinus3, nMinus1;
    
	gmp_randstate_t rnd_state;
	int found=0, finished=0;
	double start, end;
    
	mpz_init(n); mpz_init(xSquared);
    mpz_set(n, nisse);
    
	//a belongs to [1, n-2]
	//u, v <-> s belongs to [0, n-1]
	mpz_init_set(nMinus3, n);
	mpz_sub_ui(nMinus3, nMinus3, 4);
    
	mpz_init_set(nMinus1, n);
	mpz_sub_ui(nMinus1, nMinus1, 1);
    
    //[1 choose seeds]
	gmp_randinit_default(rnd_state);
	mpz_init(a); mpz_init(s);
    
	//mpz_urandomm(a, rnd_state, nMinus3);
	//mpz_add_ui(a, a, 1);
    
	//a is set to 1, comment this if you want it random
	mpz_set_ui(a, 1);
	mpz_urandomm(s, rnd_state, nMinus1);
    
	mpz_init_set(u, s);
	mpz_init_set(v, s);
    
	mpz_init(uMinusV); mpz_init(g); mpz_init_set_ui(gMul, 1);
    
	//Pollard's rho cannot tell if a number is prime, test before getting into an infinite loop
	if(mpz_probab_prime_p(n, 5) > 0)
	{
		printf("%s is prime\n", mpz_get_str(NULL, 10, n));
        mpz_set(x, n);
        mpz_set_ui(y, 1);
        return;
        //		exit(0);
	}
    
	unsigned long steps = 0;
    //	start = my_gettimeofday();
    
	while(!finished)
	{
        //[Factor search]
        while(!found)
        {
            f(u);
            f(v);
            f(v);
            
            mpz_set(uMinusV, u);
            mpz_sub(uMinusV, uMinusV, v);
            mpz_abs(uMinusV, uMinusV);
            
            //We don't calculate gcd everytime, we do 100 multiplications and use the result to
            //extract the gcd, since it must be among the product
            mpz_mul(gMul, gMul, uMinusV);
            //            if (Miller(gMul, 5)) {
            uint64_t c = 2;
            mpz_class gmultemp(gMul);
            if (Miller(gmultemp, 5) || mpz_cmp_ui(gMul, 0) == 0) {
                
                mpz_t m; mpz_init(m);
                while (true) {
                    mpz_mod_ui(m, n, c);
                    //                        printis(m); cout << " " << endl;
                    if (mpz_cmp_ui(m, 0) == 0) {
                        //                        cout << " WOO: " << c;
                        //                        cout << " " << endl;
                        mpz_set_ui(x, c);
                        mpz_t temp; mpz_init(temp);
                        mpz_div_ui(temp, n, c);
                        mpz_set(y, temp);
                        return;
                    } else {
                        c++;
                    }
                }
            }
            //            }
            if(steps%100 == 0)
            {
                //                printis(gMul);
                gcd(g, gMul, nisse);
                
                if(mpz_cmp_ui(g, 1) != 0)
                {
                    found = 1;
                }
                
                mpz_set_ui(gMul, 1);
            }
            
            steps++;
        }
        
        printf("Testing GCD: %s\n", mpz_get_str(NULL, 10, g));
        
        //[Bad seed]
        if(mpz_cmp(g, n) != 0)
            finished = 1;
	}
    
    //	end = my_gettimeofday();
    
	printf("Found divisor g = %s in %lu steps [%.3f s]\n", mpz_get_str(NULL, 0, g), steps,
           end - start);
    
    mpz_set(x, g);
    mpz_div(y, x, g);
    //    mpz_set(x, g);
    
}

int main(int argc, const char * argv[])
{
    
    
    vector<mpz_class> PRIMES;
//    allocate_primes();
    primenumbers primeobject;
    vector<int> & primes = primeobject.get_primes();
    for (int i = 0; i < primes.size(); i++) {
        mpz_class tmp(primes[i]);
        PRIMES.push_back(tmp);
    }
    primeobject.clear();
//    deallocate_primes();
    mpz_t B;
    mpz_init(B);
    
    //  mpz_init(B);
    //    mpz_set_str(N, "9207215733000000000000000000000000000000000000000000000000000000000000", 10);
    //    mpz_set_str(N, "92072157330000000000000000000000000000000", 10);
    //    mpz_set_str(N, "8741261238172833231", 10);
    //    mpz_set_str(N, "96573982938235691", 10);
    
    //    mpz_class N("96573982938235692");
    mpz_class N("9207215733000000000000000000000000000000000000000000000000000000000001");
    
    
    
    /** START OF ALGORITHM **/
    
    vector<mpz_class> result; // store resulting factors
    queue<mpz_class> cand_q; // store candidate factors
    cand_q.push(N);
fromTheTop:
    while (!cand_q.empty()) {
        
        mpz_class n;
        
        // pop element from candidate queue
        n = cand_q.front();
        cand_q.pop();
        mpz_t n_t;
        mpz_init(n_t);
        // convert mpz_class to mpz_t
        string str = n.get_str(10);
        char* p = new char[str.length()+1];
        memcpy(p, str.c_str(), str.length()+1);
        mpz_set_str(n_t, p, 10);
        
        
        
        cout << "------------------------------------" << endl;
        cout << "N: " << n;
        
        cout  << endl;
        
        if (n == 1)
            continue;
        
        //        if (Miller(n, 10000)) { // vadå klarar inte Miller av att säga att 5 är prim?
        if (mpz_probab_prime_p(n_t, 25)) {
            cout << "Storing (miller says prime): " << n;
            cout << endl;
            result.push_back(n);
        } else {
            
            quadS qs(n_t);
            qs.generateBaseLimit(B, n_t);
            int retries = 0;
            
            for (int i = 0; i < PRIMES.size(); ++i) {
                if (n/PRIMES[i] * PRIMES[i] == n) {
                    cout << "Found factorization: " << n/PRIMES[i] << " * " << PRIMES[i] << endl;
                    cand_q.push(n/PRIMES[i]);
                    cand_q.push(PRIMES[i]);
                    goto fromTheTop;
                }
            }
            
        retry:
            
            if (retries > 5 || n < 10048576) {
                cout <<"*************Pollard*************" << endl;
                
                mpz_t x_p, y_p;
                mpz_init(x_p); mpz_init(y_p);
                runPollard(n_t, x_p, y_p);
                cout << "pollard answers: ";
                printis(x_p);
                cout << " ";
                printis(y_p);
                mpz_class x_c(x_p); mpz_class y_c(y_p);
                if (x_c != n && x_c != 1) {
                    cand_q.push(x_c);
                    cand_q.push(y_c);
                }
                continue;
                
            }
            
            cout << "*************QS*************" << endl;
            cout << "Base: ";
            qs.print(B);
            
            int64_t t = qs.sieve_primes_up_to((int64_t) qs.toInt(B));
            cout << "Primes found: " << t << endl;
            int64_t w = qs.fill_primes_with_quadratic_residue(n_t);
            cout << "Quadtratic primes found: " <<  w << endl;

            qs.generateModularRoots();
            
            qs.matrix_init(w, w);
            qs.set_matrix_row(w);
            
            cout << "Start sieving" << endl;
            qs.sieve();
            cout << "Start gauss" << endl;
            qs.gauss_elimination();
            cout << "Start factorization" << endl;
            qs.factor();
            
            
            mpz_t x, y;
            mpz_init(x); mpz_init(y);
            qs.fetch_answers(x, y);
            
            cout << "(x,y): ";
            printis(x);
            cout << ", ";
            printis(y);
            cout << endl;
            
            mpz_class tmp_x(x);
            mpz_class tmp_y(y);
            
            if (tmp_x != n && tmp_x != 1) {
                cand_q.push(tmp_x);
                cand_q.push(tmp_y);
            } else {
                mpz_add(B, B, B);
                cout << "Retrying with b: ";
                printis(B);
                cout << endl;
                retries++;
                goto retry;
            }
        }
        
    }
    
    cout << endl;
    
    cout << N << " = ";
    for (int i = 0; i < result.size(); ++i) {
        cout << result[i];
        if (i != result.size()-1)
            cout << " * ";
    }
    //    }
    cout << endl;
    
    /**END**/
    
    return 0;
}

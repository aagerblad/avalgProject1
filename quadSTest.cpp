#include <iostream>
#include "pollardRho.cpp"
#include "quadS.cpp"
//#include <gmpxx.h>
using namespace std;

int g(int a, int n) {
    return (a*a + 1) % n;
}

int gcd(int a, int b)  {
    int c;
    while ( a != 0) {
        c = a;
        a = b%a;
        b = c;
    }
    return b;
}

int main(int argc, const char * argv[])
{
    
    //    qs hej("1232");
    mpz_t N, B;
    mpz_init(N);
    mpz_init(B);
    mpz_set_str(N, "87463", 10);
    
    
    quadS qs(N);
    qs.generateFactorBase(B, N);
    
    cout << "N: ";
    qs.print(N);
    cout << "Base: ";
    qs.print(B);
    
    int64_t t = qs.sieve_primes_up_to((int64_t) qs.toInt(B));
    cout << "Primes found: " << t << endl;
    int64_t w = qs.fill_primes_with_quadratic_residue(N);
    cout << "Quadtratic primes found: " <<  w << endl;
    qs.generateModularRoots();
    qs.printModularRoots();
    
    
    
    
    return 0;
}

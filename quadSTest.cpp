#include <iostream>
#include "quadS.cpp"
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
    mpz_t N, B, temp_matrix_row;
    mpz_init(N);
    mpz_init(B);
//    mpz_set_str(N, "9207215733000000000000000000000000000000000000000000000000000000000000", 10);
    mpz_set_str(N, "7", 10);
    
    
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
//    qs.printModularRoots();
    
    qs.matrix_init(w, w);
    mpz_init2(temp_matrix_row, w);
    
    cout << "Start sieving" << endl;
    qs.sieve();
    cout << "Start gauss" << endl;
    qs.gauss_elimination();
    cout << "Start factorization" << endl;
    qs.factor();
    
    
    
    return 0;
}

/* 
 * C++ Program to Implement Miller Rabin Primality Test
 */
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <gmpxx.h>
//#define ll long long
using namespace std;

void print(mpz_class a) {
//    char * s;
//    s = mpz_get_str(NULL, 10, a);
    
    cout <<  a << endl;
    
}

/*
 * calculates (a * b) % c taking into account that a * b might overflow
 */
mpz_class mulmod(mpz_class a, mpz_class b, mpz_class mod)
{
    mpz_class x, y;
//    mpz_init(x);
//    mpz_init(y);
//    mpz_mod(y, a, mod);
    y = a % mod;
//    ll x = 0,y = a % mod;
    
    while (b > 0)
    {
        
        mpz_class tmp;
//        mpz_init(tmp);
//        mpz_mod_ui(tmp, b, 2);
        tmp = b % 2;
        if (tmp == 1)
        {
//            mpz_add(x, x, y);
//            mpz_mod(x, x, mod);
//            x += y;
//            x %= mod;
            x = (x + y) % mod;
        }
//        mpz_mul_ui(y, y, 2);
//        mpz_mod(y, y, mod);
        y = (y * 2) % mod;
//        mpz_div_ui(b, b, 2);
        b /= 2;
    }
//    return x % mod;
//    mpz_mod(res, x, mod);
    return x % mod;
}
/* 
 * modular exponentiation
 */
mpz_class modulo(mpz_class base, mpz_class e, mpz_class mod)
{
    mpz_class exponent;
//    mpz_init(exponent);
//    mpz_set(exponent, e);
    exponent = e;
//    ll x = 1;
    mpz_class x, y;
//    mpz_init(x);
//    mpz_init(y);
//    mpz_set_ui(x, 1);
//    mpz_set(y, base);
    x = 1;
    y = base;
//    ll y = base;
    while (exponent > 0)
    {
        mpz_class tmp;
//        mpz_init(tmp);
//        mpz_mod_ui(tmp, exponent, 2);
        tmp = exponent % 2;
        if (tmp == 1) {
//            mpz_mul(x, x, y);
//            mpz_mod(x, x, mod);
            x = (x * y) % mod;
        }
//        mpz_mul(y, y, y);
//        mpz_mod(y, y, mod);
        y = (y * y) % mod;
//        mpz_div_ui(exponent, exponent, 2);
        exponent = exponent / 2;
    }
//    mpz_mod(res, x, mod);
    return x % mod;
}
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(mpz_class p,int iteration)
{
    
    if (p == sqrt(p) * sqrt(p)) {
//            if (sqrt(p)==floor(sqrt(p))) {
        cout << "Hej tobbe: " << sqrt(p) << endl;
        return false;
    }
    
    if (p < 2)
    {
        return false;
    }
    
    if (p != 2 && p % 2 == 0)
    {
        return false;
    }
    
    if (p == 2 || p == 3)
        return true;
    
    mpz_class s = p-1;
    while (s % 2 == 0)
    {
        s /= 2;
    }
    
    for (int i = 0; i < iteration; i++)
    {
        mpz_class a;
        
        gmp_randstate_t r;
        gmp_randinit_default(r);
        gmp_randseed_ui(r, rand());
        
        gmp_randclass rr(gmp_randinit_default);
        
        rr.seed(time(NULL));
        a = rr.get_f();
        a = (a % (p-1)) + 1;
        mpz_class temp = s;
        mpz_class mod = modulo(a, temp, p);
        
        while (temp != p-1 && mod != 1 && mod != p-1)
        {
            mod = mulmod(mod, mod, p);
            temp *= 2;
        }
        
        if (mod != p-1 && temp % 2 == 0)
        {
            return false;
        }
    }
    return true;
}
// //Main
// int main()
// {
//     int iteration = 5;
//     ll num;
//     cout<<"Enter integer to test primality: ";
//     cin>>num;
//     if (Miller(num, iteration))
//         cout<<num<<" is prime"<<endl;
//     else
//         cout<<num<<" is not prime"<<endl;
//     return 0;
// }
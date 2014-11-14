/* 
 * C++ Program to Implement Miller Rabin Primality Test
 */
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <gmpxx.h>
//#define ll long long
using namespace std;
 
/* 
 * calculates (a * b) % c taking into account that a * b might overflow 
 */
void mulmod(mpz_t res, mpz_t a, mpz_t b, mpz_t mod)
{
    mpz_t x, y;
    mpz_init(x);
    mpz_init(y);
    mpz_mod(y, a, mod);
//    ll x = 0,y = a % mod;
    
    while (b > 0)
    {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mod_ui(tmp, b, 2);
        
        if (mpz_cmp_ui(tmp, 1) == 0)
        {
            mpz_add(x, x, y);
            mpz_mod(x, x, mod);
//            x = (x + y) % mod;
        }
        mpz_mul_ui(y, y, 2);
        mpz_mod(y, y, mod);
//        y = (y * 2) % mod;
        mpz_div_ui(b, b, 2);
//        b /= 2;
    }
//    return x % mod;
    mpz_mod(res, x, mod);
}
/* 
 * modular exponentiation
 */
void modulo(mpz_t res, mpz_t base, mpz_t exponent, mpz_t mod)
{
//    ll x = 1;
    mpz_t x, y;
    mpz_init(x);
    mpz_init(y);
    mpz_set_ui(x, 1);
    mpz_set(y, base);
//    ll y = base;
    while (exponent > 0)
    {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mod_ui(tmp, exponent, 2);
        if (mpz_cmp_ui(tmp, 1) == 0) {
            mpz_mul(x, x, y);
            mpz_mod(x, x, mod);
//            x = (x * y) % mod;
        }
        mpz_mul(y, y, y);
        mpz_mod(y, y, mod);
//        y = (y * y) % mod;
        mpz_div_ui(exponent, exponent, 2);
//        exponent = exponent / 2;
    }
    mpz_mod(res, x, mod);
//    return x % mod;
}
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(mpz_t p,int iteration)
{
    
    if (mpz_cmp_ui(p, 2) < 0)
    {
        return false;
    }
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod_ui(tmp, p, 2);
    if (mpz_cmp_ui(p, 2) != 0 && mpz_cmp_ui(tmp, 0) == 0)
    {
        return false;
    }
    mpz_t s;
    mpz_init(s);
    mpz_sub_ui(s, p, 1);
//    ll s = p - 1;
    mpz_mod_ui(tmp, s, 2);
    while (mpz_cmp_ui(tmp, 0) == 0)
    {
        mpz_div_ui(s, s, 2);
//        s /= 2;
    }
    for (int i = 0; i < iteration; i++)
    {
//        ll a = rand() % (p - 1) + 1, temp = s;
        mpz_t a;
        mpz_init(a);
        gmp_randstate_t r;
        gmp_randinit_default(r);
        mpz_urandomm(a, r, p);
        mpz_set(tmp, s);
        mpz_t mod;
        mpz_init(mod);
        modulo(mod, a, tmp, p);
//        ll mod = modulo(a, temp, p);
        
        mpz_t p_1; mpz_init(p_1);
        mpz_sub_ui(p, p, 1);
        bool arg1 = mpz_cmp(tmp, p_1) != 0;
        bool arg2 = mpz_cmp_ui(mod, 1) != 0;
        bool arg3 = mpz_cmp(mod, p_1) != 0;
        while (arg1 && arg2 && arg3)
        {
            mulmod(mod, mod, mod, p);
//            mod = mulmod(mod, mod, p);
            mpz_mul_ui(tmp, tmp, 2);
//            temp *= 2;
        }
        arg1 = mpz_cmp(mod, p_1) != 0;
        mpz_t tmpMod2; mpz_init(tmpMod2);
        mpz_mod_ui(tmpMod2, tmp, 2);
        arg2 = mpz_cmp_ui(tmpMod2, 0) == 0;
        if (arg1 && arg2)
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
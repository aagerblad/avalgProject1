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
void mulmod(mpz_class res, mpz_class a, mpz_class b, mpz_class mod)
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
    res = x % mod;
}
/* 
 * modular exponentiation
 */
void modulo(mpz_class res, mpz_class base, mpz_class e, mpz_class mod)
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
    res = x % mod;
}
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(mpz_class p,int iteration)
{
    
    if (p < 2)
    {
        return false;
    }
    
//    if (p == 3) {
//        return false;
//    }
    mpz_class tmp;
//    mpz_init(tmp);
//    mpz_mod_ui(tmp, p, 2);
    tmp = p % 2;
    if (p != 2 && tmp == 0)
    {
        return false;
    }
    mpz_class s;
//    mpz_init(s);
//    mpz_sub_ui(s, p, 1);
    s = p - 1;
//    mpz_mod_ui(tmp, s, 2);
    tmp = s % 2;
    while (tmp == 0)
    {
//        mpz_div_ui(s, s, 2);
//        mpz_mod_ui(tmp, s, 2);
        s /= 2;
        tmp = s % 2;
    }
    for (int i = 0; i < iteration; i++)
    {
//        ll a = rand() % (p - 1) + 1, temp = s;
        mpz_class a;
//        mpz_init(a);
        gmp_randstate_t r;
        gmp_randinit_default(r);
        unsigned long int seed;
        seed = rand();
        gmp_randseed_ui(r, seed);
        
//        mpz_class ran;
        gmp_randclass rr(gmp_randinit_default);
        rr.seed(time(NULL));
        a = rr.get_f();
//        ran =rr.get_z_bits(125);
//        long int random=ran.get_ui();
        
//        mpz_urandomm(a, r, p);
//        mpz_set(tmp, s);
        tmp = s;
        mpz_class mod;
//        mpz_init(mod);
        modulo(mod, a, tmp, p);
//        print(tmp);
//        ll mod = modulo(a, temp, p);
        mpz_class p_1; // mpz_init(p_1);
//        mpz_sub_ui(p_1, p, 1);
        p_1 = p - 1;
//        bool arg1 = mpz_cmp(tmp, p_1) != 0;
//        bool arg2 = mpz_cmp_ui(mod, 1) != 0;
//        bool arg3 = mpz_cmp(mod, p_1) != 0;
        while (tmp != p_1 && mod != 1 && mod != p_1)
        {
            mulmod(mod, mod, mod, p);
//            mod = mulmod(mod, mod, p);
//            mpz_mul_ui(tmp, tmp, 2);
            tmp *= 2;
//            cout << "mod: ";
//            print(mod);
//            cout << "tmp: ";
//            print(tmp);
//            arg1 = mpz_cmp(tmp, p_1) != 0;
//            arg2 = mpz_cmp_ui(mod, 1) != 0;
//            arg3 = mpz_cmp(mod, p_1) != 0;
        }
//        arg1 = mpz_cmp(mod, p_1) != 0;
        mpz_class tmpMod2; //mpz_init(tmpMod2);
        //mpz_mod_ui(tmpMod2, tmp, 2);
        tmpMod2 = tmp % 2;
//        arg2 = mpz_cmp_ui(tmpMod2, 0) == 0;
        if (mod != p_1 && tmpMod2 == 0)
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
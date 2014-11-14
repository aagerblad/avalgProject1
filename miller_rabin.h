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
void mulmod(mpz_t res, mpz_t a, mpz_t b, mpz_t mod);
/* 
 * modular exponentiation
 */
void modulo(mpz_t res, mpz_t base, mpz_t exponent, mpz_t mod);
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(mpz_t p, int iteration);
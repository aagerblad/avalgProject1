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
void mulmod(mpz_class res, mpz_class a, mpz_class b, mpz_class mod);
/* 
 * modular exponentiation
 */
void modulo(mpz_class res, mpz_class base, mpz_class exponent, mpz_class mod);
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(mpz_class p, int iteration);
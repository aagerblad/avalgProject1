/* 
 * C++ Program to Implement Miller Rabin Primality Test
 */
#include <iostream>
#include <cstring>
#include <cstdlib>
#define ll long long
using namespace std;
 
/* 
 * calculates (a * b) % c taking into account that a * b might overflow 
 */
ll mulmod(ll a, ll b, ll mod);
/* 
 * modular exponentiation
 */
ll modulo(ll base, ll exponent, ll mod);
 
/*
 * Miller-Rabin primality test, iteration signifies the accuracy
 */
bool Miller(ll p,int iteration);
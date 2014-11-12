#include <iostream>
#include "pollardRho.cpp"
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
    
 //    int n = 10000;
	// int x_fixed = 4711;
	// int cycle_size = 2;
	// int x = x_fixed;
	// int h = 1;
    
	// while (h == 1) {
	// 	int count = 1;
        
	// 	while (count <= cycle_size && h == 1) {
	// 		x = g(x, n);
	// 		count = count + 1;
	// 		h = gcd(x - x_fixed, n);
	// 	}
        
	// 	if (h != 1)
	// 		break;
        
	// 	cycle_size = 2 * cycle_size;
	// 	x_fixed = x;
	// }
	// cout << "\nThe factor is  " << h;
    
	pollardRho p(10050001);

	std::vector<long long> v = p.getAllFactors();

	for (auto it = v.begin(); it != v.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;

	// int iteration = 5;
	// if (Miller(17, 5)) {
	// 	cout << 17 << " is prime" << endl;
 // 	} else {
 // 		cout << 17 << " is not prime" << endl;
 // 	}

	// cout << p.gcd(4712, 1244) << endl;

	// cout << p.g(5) << endl;
   
   return 0;
}

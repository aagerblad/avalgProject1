#include <iostream>
#include "pollardRho.cpp"
using namespace std;

int main(int argc, const char * argv[])
{
	pollardRho p(10403);

	// cout << p.calculateFactor() << endl;

	std::vector<long long> v = p.getAllFactors();

	for (auto it = v.begin(); it != v.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;

	// cout << p.gcd(4712, 1244) << endl;

	// cout << p.g(5) << endl;
   
   return 0;
}

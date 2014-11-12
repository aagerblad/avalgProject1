#include <iostream>
#include <vector>
using namespace std;

class pollardRho
{
private:
	long long startValue;
	

public:
	pollardRho(long long n_n) {
		startValue = n_n;
	}
	~pollardRho() {}

	std::vector<long long> getAllFactors() {
		std::vector<long long> v;
		long long n = startValue;
		long long newStart = n;
		long long x_0 = 2;
		while (n != 1) {
			long long factor = calculateFactor(n, x_0);
			if (factor < 1) {
				// v.clear();	
				// n = newStart;
				++x_0;
				// cout << x_0 << endl;
			} else {
				v.push_back(factor);
				n = n/factor;
				// newStart = n;
				x_0 = 2;
			}
		}
		cout << x_0 << endl;
		return v;
	}

	long long calculateFactor(long long n, long long x_0) {
		long long rValue = 1;
		// long long x_0 = 4;
		long long x_i = x_0;
		long long cycle_size = 2;

		// x_i = g(x_0);
		// rValue = gcd(x_i - x_0, n);
		
		while (rValue == 1) {
			
			long long count = 1;
			while (count <= cycle_size && rValue == 1) {
				++count;
				x_i = g(n, x_i);
				rValue = gcd(x_i - x_0, n);
			}
	
			if (rValue != 1)
				break;
	
			cycle_size *= 2;
			x_0 = x_i;

		}
	
		return rValue;
	}

	long long g(long long n, long long a) const {
		return (a*a + 1) % n;
	}

	long long gcd(long long a, long long b) const {
		long long c;
		while ( a != 0) {
			c = a;
			a = b%a;
			b = c;
		}
		return b;
	}
};
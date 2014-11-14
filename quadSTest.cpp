#include <iostream>
#include <queue>
#include <gmpxx.h>
#include "quadS.cpp"
#include "miller_rabin.h"
using namespace std;

int hej = 0;

void printis(mpz_t a) {
    char * s;
    s = mpz_get_str(NULL, 10, a);
    cout <<  s;
    
}

class mpzqueue {
    
private:
    mpz_class * elements;
    int totalelements;
//    queue<node> queue;
    
public:
    
    mpzqueue() {
        elements = new mpz_class[200];//(mpz_t *) calloc(200, sizeof(mpz_t));
        totalelements = 0;
    }
    void push(mpz_class d) {
        hej++;
//        mpz_t e; mpz_init(e);
//        mpz_set(e, d);
//        mpz_set(elements[totalelements], e);
        elements[totalelements] = d;
        totalelements++;
    }
    
    void pop(mpz_class ret) {
        if (totalelements<= 0) {
            return;
        }
        hej--;

//        mpz_set(ret, elements[totalelements-1]);
//        mpz_set_ui(elements[totalelements-1], 0);
        ret = elements[totalelements-1];
        elements[totalelements-1] = 0;
        totalelements--;
        
    }
    
    bool isEmpty() {
        if (totalelements <= 0)
            return true;
        return false;

    }
    
    void print() {
        cout << "size: " << totalelements << endl;
        for (int i = 0; i < totalelements; i++) {
            cout << elements[i] << " ";
//            if (mpz_cmp_ui(elements[i], 0) != 0) {
//                printis(elements[i]);
//                cout << " ";
//            }
        }
        cout << endl;
        cout << "HEJ: " << hej << endl;
    }
};

void gcd(mpz_t & ret, mpz_t a, mpz_t b) {
    mpz_t c;
    mpz_init(c);
    mpz_t zero;
    mpz_init(zero);
    while (mpz_cmp(a, zero) != 0) {
        mpz_set(c, a);
        mpz_mod(a, b, a);
        mpz_set(b, c);
    }
    mpz_set(ret, b);
}


mpz_t xSquared, a, n;

//rop <- (rop² + a) mod n
void f(mpz_t rop)
{
	mpz_pow_ui(xSquared, rop, 2);
	mpz_add(xSquared, xSquared, a);
	mpz_mod(rop, xSquared, n);
}

void runPollard(mpz_t nisse, mpz_t x, mpz_t y)
{
//	mpz_t g, gMul, s, u, v, uMinusV;
//	mpz_t nMinus3, nMinus1;
//    
//	gmp_randstate_t rnd_state;
//	int found=0, finished=0;
//	double start, end;
//    
//	mpz_init(n); mpz_init(xSquared);
//    mpz_set(n, nisse);
//    
//	//a belongs to [1, n-2]
//	//u, v <-> s belongs to [0, n-1]
//	mpz_init_set(nMinus3, n);
//	mpz_sub_ui(nMinus3, nMinus3, 4);
//    
//	mpz_init_set(nMinus1, n);
//	mpz_sub_ui(nMinus1, nMinus1, 1);
//    
//    //[1 choose seeds]
//	gmp_randinit_default(rnd_state);
//	mpz_init(a); mpz_init(s);
//    
//	//mpz_urandomm(a, rnd_state, nMinus3);
//	//mpz_add_ui(a, a, 1);
//    
//	//a is set to 1, comment this if you want it random
//	mpz_set_ui(a, 1);
//	mpz_urandomm(s, rnd_state, nMinus1);
//    
//	mpz_init_set(u, s);
//	mpz_init_set(v, s);
//    
//	mpz_init(uMinusV); mpz_init(g); mpz_init_set_ui(gMul, 1);
//    
//	//Pollard's rho cannot tell if a number is prime, test before getting into an infinite loop
//	if(mpz_probab_prime_p(n, 5) > 0)
//	{
//		printf("%s is prime\n", mpz_get_str(NULL, 10, n));
//        mpz_set(x, n);
//        mpz_set_ui(y, 1);
//        return;
////		exit(0);
//	}
//    
//	unsigned long steps = 0;
//    //	start = my_gettimeofday();
//    
//	while(!finished)
//	{
//        //[Factor search]
//        while(!found)
//        {
//            f(u);
//            f(v);
//            f(v);
//            
//            mpz_set(uMinusV, u);
//            mpz_sub(uMinusV, uMinusV, v);
//            mpz_abs(uMinusV, uMinusV);
//            
//            //We don't calculate gcd everytime, we do 100 multiplications and use the result to
//            //extract the gcd, since it must be among the product
//            mpz_mul(gMul, gMul, uMinusV);
////            if (Miller(gMul, 5)) {
//            uint64_t c = 2;
//                if (Miller(gMul, 5) || mpz_cmp_ui(gMul, 0) == 0) {
//                    
//                    mpz_t m; mpz_init(m);
//                    while (true) {
////                        cout << "vafa"
////                        unsigned long d = mpz_get_ui(nisse);
////                        cout << "c: " << c << " d: "  << d << endl;
//                        
//                        mpz_mod_ui(m, n, c);
////                        printis(m); cout << " " << endl;
//                    if (mpz_cmp_ui(m, 0) == 0) {
////                        cout << " WOO: " << c;
////                        cout << " " << endl;
//                        mpz_set_ui(x, c);
//                        mpz_t temp; mpz_init(temp);
//                        mpz_div_ui(temp, n, c);
//                        mpz_set(y, temp);
//                        return;
//                    } else {
//                     c++;
//                    }
//                    }
//                }
////            }
//            if(steps%100 == 0)
//            {
////                printis(gMul);
//                gcd(g, gMul, nisse);
//                
//                if(mpz_cmp_ui(g, 1) != 0)
//                {
//                    found = 1;
//                }
//                
//                mpz_set_ui(gMul, 1);
//            }
//            
//            steps++;
//        }
//        
//        printf("Testing GCD: %s\n", mpz_get_str(NULL, 10, g));
//        
//        //[Bad seed]
//        if(mpz_cmp(g, n) != 0)
//            finished = 1;
//	}
//    
//    //	end = my_gettimeofday();
//    
//	printf("Found divisor g = %s in %lu steps [%.3f s]\n", mpz_get_str(NULL, 0, g), steps,
//           end - start);
//    
//    mpz_set(x, g);
//    mpz_div(y, x, g);
////    mpz_set(x, g);
    
}

int main(int argc, const char * argv[])
{
    
    
//    mpzqueue cq;
//    
//    mpz_t hej, hej2, hej3, hej4;
//    mpz_init(hej); mpz_init(hej2); mpz_init(hej3); mpz_init(hej4);
//    mpz_set_ui(hej, 1);
//        mpz_set_ui(hej2, 2);
//        mpz_set_ui(hej3, 3);
//        mpz_set_ui(hej4, 4);
//    
//    cout << "pushar" << endl;
//    cq.push(hej);
//    cq.print();
//    cout << "pushar" << endl;
//    cq.push(hej2);
//    cq.print();
//    cout << "poppar" << endl;
//    cq.pop(hej3);
//    cq.print();
//    cout << endl;
//    cout << "hej3: ";
//    printis(hej3);
//    cout << endl;
//    cq.push(hej4);
//    cq.print();
//    cout << "SLUT PÅ TEST"  << endl;

    
    //    qs hej("1232");
    mpz_t B, temp_matrix_row;
//    mpz_init(N);
    mpz_init(B);
    //    mpz_set_str(N, "9207215733000000000000000000000000000000000000000000000000000000000000", 10);
    //    mpz_set_str(N, "92072157330000000000000000000000000000000", 10);
//    mpz_set_str(N, "8741261238172833231", 10);
//    mpz_set_str(N, "96573982938235691", 10);
    
    mpz_class N("96573982938235691");
    
    /**START**/

//    mpz_t * result = new mpz_t[100];
    vector<mpz_class> result;

//    mpzqueue cand_queue;
//    cand_queue.push(N);
//    cout << "first: ";
//    cand_queue.print();
//    mpz_set(n, N);
    int added_results = 0;
    mpzqueue cand_q;
    cand_q.push(N);
//    mpz_t * cands = new mpz_t[300];
//    mpz_set(cands[0], N);
//    int candsize = 1;
//    int current_cand = 0;
    while (!cand_q.isEmpty()) {

        mpz_class n;
//        cout << "start of loop: ";
//        cand_queue.print();
        
//        mpz_set(n, cand_q.);
          cand_q.pop(n);
//        mpz_init(n);
        
//        cand_queue.pop(n);
        
//        cout << "after first pop: ";
//        cand_queue.print();
//        cout << "n efter första popen: ";
//        printis(n);
        
//        mpz_t b; mpz_init(b);
//        mpz_set_ui(b, 97);
//        if(Miller(b, 5)) {
//            // TODO varför tror den att detta tal är prim?
//            cout << "NÄMEN VAFAN" << endl;
//        }
    
        cout  << endl;
//        printis(n);
        if (Miller(n, 5)) { // vadå klarar inte Miller av att säga att 5 är prim?
//            cout << "halloj" << endl;
            
            
            cout << "Storing: " << n;
//            printis(n);
            cout << endl;
//            --current_cand;
            //            results.push_back(n);
//            cand_queue.pop(result[added_cands]);
//            mpz_set(result[added_results],n);
//            results[added_results], n);
            result.push_back(n);
//            vector<mpz_class> hej;
            
            added_results++;
//            current_cand++;
        } else {
            
            
            
            //        bool m =
            //        cout << "ÄRE PRIM? : " << m << endl;
            quadS qs(n);
            cout << "------------------------------------" << endl;
            cout << "N: ";
            qs.print(n);
            
            if (n < 0) {
                // TODO kanske köra pollard för små tal?
//                cout << "*************Pollard*************" << endl;
//                mpz_t ena, andra;
//                mpz_init(ena); mpz_init(andra);
//                runPollard(n, ena, andra);
//                cout << "två tal: ";
//                printis(ena);
//                cout << " ";
//                printis(andra);
//                cout << endl;
//                mpz_set(cands[candsize], ena);
//                candsize++;
//                mpz_set(cands[candsize], andra);
//                candsize++;
//                current_cand++;
//                cand_queue.push(ena);
//                cout << "after push ena:";
//                cand_queue.print();
//                cand_queue.push(andra);
//                cout << "afer push andra:";
//                cand_queue.print();
                
            } else {
                
            cout << "*************QS*************" << endl;
            qs.generateBaseLimit(B, n);
            
            cout << "Base: ";
            qs.print(B);
            
            int64_t t = qs.sieve_primes_up_to((int64_t) qs.toInt(B));
            cout << "Primes found: " << t << endl;
            int64_t w = qs.fill_primes_with_quadratic_residue(n);
            cout << "Quadtratic primes found: " <<  w << endl;
            //    qs.printPrimes();
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
            
            mpz_t x, y;
            mpz_init(x); mpz_init(y);
            qs.fetch_answers(x, y);

            cout << "(x,y): ";
            printis(x);
                cout << ", ";
            printis(y);
                cout << endl;
                
                current_cand++;
                
                mpz_set(cands[candsize], x);
                candsize++;
                mpz_set(cands[candsize], y);
                candsize++;
//                cand_queue.print();
//                cand_queue.push(x);
//                cout << "after next to last push:";
//                cand_queue.print();
//                cand_queue.push(y);
//                cout << "after last push:";
//                cand_queue.print();
  
            }
        }
        
        
    }
    
    cout << endl;
    cout << mpz_get_str(NULL, 10, N) << " = ";
    for (int i = 0; i < 100; i++) {
        if (mpz_cmp_ui(result[i], 0) != 0) {
            printis(result[i]);
        } else {
            break;
        }
        if (mpz_cmp_ui(result[i+1], 0) != 0) {
            cout << " * ";
        }
    }
    cout << endl;
    
    /**END**/
    
    
    
    
    return 0;
}

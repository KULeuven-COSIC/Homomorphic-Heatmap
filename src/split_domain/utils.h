#ifndef __MY_UTILS_FUNCS__
#define __MY_UTILS_FUNCS__

//#include <NTL/ZZ.h>
//#include <NTL/ZZ_p.h>
//#include <NTL/ZZX.h>
//#include <NTL/mat_ZZ_p.h>
//#include <NTL/mat_ZZ.h>
//#include <NTL/RR.h>

#include "flint/fmpz_polyxx.h"

#include "fft.h"

#include <vector>
#include <iostream> 

using namespace std;
using namespace flint;
//using namespace NTL;
//

bool rand_bit();


// returns inverse of a mod n
long inverse(long a, long n);

// Let P = p[0] * p[1] * ... * p[n-1] and Pi := P / p[i].
// At the end of function, we have invP[o] = Pi^-1 mod p[i].
void compute_inverses_crt(std::vector<int>& invP, const std::vector<int>& p);

// random degree N-1 polynomial with coefficients in [[0, 2^logq - 1]]
fmpz_polyxx random_poly(long N, long logq);

//void random_poly32(Poly32& a, int N, int logn);

fmpzxx symmetric_mod(const fmpzxx& a, const fmpzxx& n); // returns a % n using the set ]-c/2, c/2] as Z/nZ

void reduce_poly_mod(fmpz_polyxx& poly, fmpzxx mod, bool centered = true);

// Receives a polynomial of degree smaller than 2N and reduces it mod x^N+1
void reduce_2N_to_N(fmpz_polyxx& poly, long N);
/**
 * Return  round(a/n), that is, interpret a/n as an rational number
 * then return the closest integer to it.
 */
fmpzxx rounded_division(const fmpzxx& a, const fmpzxx& n);

fmpz_polyxx rounded_division(const fmpz_polyxx& a, const fmpzxx& n);

fmpzxx max_norm(const fmpz_polyxx& poly);

template <typename T1, typename T2>
vector<T1>& operator+= (vector<T1>& vec, const vector<T2>& u);

/**    Set words_vec to the vector (x0, x1, ..., x_{l-1})
 * representing the decomposition of x in base b.
 */
void decompose_ZZ(std::vector<fmpzxx>& words_vec, const fmpzxx& x, long l, const fmpzxx& b);

/**
 *	Assumes that deg(f) < N and that all the coefficients of f are smaller than
 * b^l in absolute value.
 * 	Sets words_vec to the vector (y0, y1, ..., y_{l-1})
 * representing the decomposition of f in base b, that is, each yi is a
 * polynomial with degree smaller than N and coefficients in ]-b, b[, and
 * the sum y0*b^0 + y1*b^1 + ... + y_{l-1}*b^{l-1} equals f
 */
void decompose_ZZX(std::vector<fmpz_polyxx>& words_vec, const fmpz_polyxx& f, long N, long l, const fmpzxx& b);


template <typename T>
std::ostream& operator<< (std::ostream &out, const vector<T> & u) {
	if (0 == u.size())
		return out << "[ ]";
	cout << "[";
	for (long i = 0; i < u.size()-1; i++)
		out << u[i] << ", ";
	out << u[u.size()-1] << "]";
	return out;
}



/**
 *	Returns a(X^z) mod X^N + 1
 */ 
fmpz_polyxx eval_poly(const fmpz_polyxx& a, long z, long N);

fmpz_polyxx Poly64_to_PolyFlint(const Poly64& poly);

Poly64 polyFlint_to_Poly64(const fmpz_polyxx& poly);




#endif

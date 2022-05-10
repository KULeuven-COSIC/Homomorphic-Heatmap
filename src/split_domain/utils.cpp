#include "utils.h"

#include <cassert>
using namespace std;
//using namespace NTL;
using namespace flint;


frandxx state;


bool rand_bit(){
	return (rand() % 2);
}

// returns inverse of a mod n
long inverse(long a, long n){
	unsigned long int inv;
	n_gcdinv(&inv, a%n, n);
	return inv;	
}

void compute_inverses_crt(vector<int>& invP, const vector<int>& p){
	fmpzxx P(p[0]);
	for(int i = 1; i < p.size(); i++){
		P *= p[i];
	}
	for(int i = 0; i < p.size(); i++){
		fmpzxx tmp(P / p[i]);
		tmp %= fmpzxx(p[i]);
		long Pi = tmp.to<long>();
		invP[i] = (int)inverse(Pi, (long)p[i]);
	}
}


// Receives a polynomial of degree smaller than 2N and reduces it mod x^N+1
void reduce_2N_to_N(fmpz_polyxx& poly, long N){
	for(long i = N; i < 2*N && i <= poly.degree(); i++){
		poly.coeff(i-N) -= poly.get_coeff(i);
		poly.set_coeff(i, 0);
	}
}


fmpz_polyxx random_poly(long N, long logq){
	fmpz_polyxx a;
	fmpzxx ai;
	for (long i = N-1; i >= 0; i--){
		ai = fmpzxx::randtest(state, logq-1);
		if(rand_bit())
			a.set_coeff(i, ai);
		else
			a.set_coeff(i, -ai);
	}
	return a;
}
//
//void random_poly32(Poly32& a, int N, int logn){
//	if (a.size() != N){
//		a.resize(N);
//	}
//	int32_t n_minus_1 = ((int32_t)(1 << logn)) - 1;
//	for (long i = 0; i < N; i++){
//		a[i] = rand() ; // XXX: not cryptographically secure
//		a[i] &= n_minus_1; // reduce mod n
//	}
//}
fmpzxx symmetric_mod(const fmpzxx& a, const fmpzxx& n){ // returns a % n using the set ]-c/2, c/2] as Z/nZ
	fmpzxx b(a % n);
	if(2*b >= n)
		b = b - n;
	return b;
}

void reduce_poly_mod(fmpz_polyxx& poly, fmpzxx mod, bool centered){
	fmpzxx ai;
	for (long i = 0; i <= poly.degree(); i++){
		ai = poly.get_coeff(i);
		if (centered)
			poly.set_coeff(i, symmetric_mod(ai, mod));
		else
			poly.set_coeff(i, ai % mod);
	}
}

/**
 * Return  round(a/n), that is, interpret a/n as an rational number
 * then return the closest integer to it.
 */
fmpzxx rounded_division(const fmpzxx& a, const fmpzxx& n){
	fmpzxx signal(a >= 0 ? 1 : -1);
	fmpzxx _a(a*signal);
	// interpret a = q*n + r with 0 <= r < n
	fmpzxx q(_a / n);
	fmpzxx r(_a % n);
	if (2*r < n)
		return fmpzxx(signal * q);
	else
		return fmpzxx(signal * (q + 1));
}

/** Let f = sum_{i=0}^N f_i * x^i be a polynomial of degree n.
 *  This function returns
 *      sum_{i=0}^N round(f_i / n) * x^i
 * that is, it applies rounded_division(ZZ, ZZ) to each coefficient of f.
 */
fmpz_polyxx rounded_division(const fmpz_polyxx& f, const fmpzxx& n){
	fmpz_polyxx result;
	fmpzxx f_i;
	for (long i = f.degree(); i >= 0; i--){
		f_i = f.get_coeff(i);
		result.set_coeff(i, rounded_division(f_i, n));
	}
	return result;
}


fmpzxx max_norm(const fmpz_polyxx& poly){
	fmpzxx inf_norm(0);
	for(long i = 0; i < poly.degree(); i++){
		fmpzxx abs_ai(abs(poly.get_coeff(i)));
		if (inf_norm < abs_ai)
			inf_norm = abs_ai;
	}
	return inf_norm;
}

template <typename T1, typename T2>
vector<T1>& operator+= (vector<T1>& vec, const vector<T2>& u){
//vector<fmpz_polyxx>& operator+= (vector<fmpz_polyxx>& vec, const vector<fmpz_polyxx>& u){
	for (long i = 0; i < vec.size(); i++)
		vec[i] += u[i];
	return vec;	
}


/**    Set words_vec to the vector (x0, x1, ..., x_{l-1})
 * representing the decomposition of x in base b.
 */
void decompose_ZZ(std::vector<fmpzxx>& words_vec, const fmpzxx& _x, long l, const fmpzxx& b){
	fmpzxx x(_x);
	long sign = (x < 0 ? -1 : 1);
	x *= sign;
	for (int j = 0; j < l; j++){
		words_vec[j] = sign * (x % b);
		x /= b;
	}
}

/**    Set words_vec to the vector (y0, y1, ..., y_{l-1})
 * representing the decomposition of f in base b, that is, each yi is a
 * polynomial with degree smaller than N and coefficients in ]-b, b[, and
 * the sum y0*b^0 + y1*b^1 + ... + y_{l-1}*b^{l-1} equals f
 */
void decompose_ZZX(std::vector<fmpz_polyxx>& words_vec, const fmpz_polyxx& f, long N, long l, const fmpzxx& B){
	fmpz_polyxx pow_x(1);
	fmpz_polyxx _x(0); 
	pow_x.set_coeff(0, fmpzxx(1));
	_x.set_coeff(1, fmpzxx(1));
	vector<fmpzxx> decomp_coef(l);
	fmpzxx fi;
	// decomp_coef = g^-1(f_0)
	fi = f.get_coeff(0);
	decompose_ZZ(decomp_coef, fi, l, B);
	for (long j = 0; j < l; j++)
		words_vec[j] = pow_x * decomp_coef[j];

	for(long i = 1; i < N; i++){
		pow_x *= _x; // pow_x = x^i 
		fi = f.get_coeff(i);
		if (fi != 0){
			decompose_ZZ(decomp_coef, fi, l, B);
			// After this loop: words_vec += x^i * g^-1(f_i)
			for (long j = 0; j < l; j++){
				words_vec[j] += pow_x * decomp_coef[j];
			}
		}
	}
}


// Returns a(X^z) mod X^N + 1
fmpz_polyxx eval_poly(const fmpz_polyxx& a, long z, long N){
	fmpz_polyxx b;
	for (long i = 0; i <= a.degree(); i++){
		long r = i*z % N;
		long q = (i*z - r) / N; // i*z = q*N + r
		if (0 == q % 2){
			b.set_coeff(r, a.get_coeff(i));
		}else{
			b.set_coeff(r, -a.get_coeff(i));
		}
	}
	return b;
}

fmpz_polyxx Poly64_to_PolyFlint(const Poly64& poly){
	fmpz_polyxx a;
	for(int i = 0; i < poly.size(); i++){
		a.set_coeff(i, poly[i]);
		assert(poly[i] == a.get_coeff(i));	

	}
	return a;
}


Poly64 polyFlint_to_Poly64(const fmpz_polyxx& poly){
	Poly64 a(poly.degree()+1);
	for(int i = 0; i < a.size(); i++){
		a[i] = poly.get_coeff(i).to<long>();
		assert(a[i] == poly.get_coeff(i));	
	}
	return a;
}



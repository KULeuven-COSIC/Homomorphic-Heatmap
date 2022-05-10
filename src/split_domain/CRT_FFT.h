#ifndef __CRT_FFT__
#define __CRT_FFT__

#include <vector>
//#include <fftw3.h>
#include "fft.h"
#include <complex>
#include <map>
#include <cassert>

#include "flint/fmpzxx.h"

#include "PolyCRT_FFT.h"
#include "utils.h"



/**
 * 		Class to perform exact multiplication of polynomials using FFTs instead 
 * 	of NTTs. To avoid the approximation erros introduced by the FFT, the 
 * 	polynomials are first decomposed modulo some primes p1, ..., pl,
 * 	then l FFTs are performed. Since each prime pi is small, the result mod pi
 * 	has an error smaller than 1/2, thus, it is erased by rounding.
 * 		Once we have the exact result mod each pi, we can obtain the expected
 * 	result by applying the CRT to each coefficient.
 *
 * 		This class is designed with the following application in mind:
 * 			- Let N be a power of two;
 * 			- Let R = ZZ[X] / <X^N + 1> be a cyclotomic (polynomial) ring;
 * 			- Let Rl be the set of polynomials of R whose coefficients have bl bits;
 * 			- Let Rr be the set of polynomials of R whose coefficients have br bits;
 * 			- Let u be a d-dimensional vector with entries u[0],...,u[d-1] in Rl;
 * 			- Let v be a d-dimensional vector with entries v[0],...,v[d-1] in Rr;
 * 		then, using this class, we can compute the inner-product u*v \in R.
 */
class CRT_FFT
{
 
	public:

		FFT_engine* fft;

		int logN;
		int N;

		int d_max;  // inner products are computed by chunks of d_max elements

		int bl; // bit length of coefficients of left operand

		int br; // bit length of coefficients of right operand

		// CRT-related values
		int l; // number of primes
		std::vector<int> p; // 20-bit primes
		std::vector<int> invP; // Let P = p1 * ... * pl, then invP[i] = (P/pi)^-1 mod pi
		flint::fmpzxx P; // p1 * ... * pl
		vector<flint::fmpzxx> Phat; // Phat[i] = P / pi
		vector<flint::fmpzxx> Phat_invP; // Phat[i] = (P / pi) * ((P/pi)^-1 mod pi) mod P

	/** 
	 *	logN is the logarithm of N in base 2
	 *	bl: bit length of coefficients of the polynomials of the left vector
	 *	br: bit length of coefficients of the polynomials of the right vector
	 *	d is the dimension of the vectors that will be multiplied
	 */
	CRT_FFT(int logN, int bl, int br, int d);

	~CRT_FFT();

	template <typename T> // T is supposed to be int32_t, int64_t, ZZ...
	std::vector<Poly32> reduce(const vector<T>& poly){
		std::vector<Poly32> residues(this->l);
		for(int i = 0; i < l; i++){
			residues[i] = Poly32(N, 0);
			for(int j = 0; j < poly.size(); j++){
				residues[i][j] = poly[j] % p[i];
			}
		}
		return residues;
	}

	// This function is the inverse of reduce, that is, it takes a vector
	// of polynomials mod p1,...,pl and return the original polynomial 
	// mod p1 * ... * pl
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> CRT(const std::vector<Poly32>& vec_polys){
		vector<T> ans(N);
		flint::fmpzxx mod(1L << 62);
		mod *= 4;
		for(int i = 0; i < N; i++){
			flint::fmpzxx ai(0);
			for(int j = 0; j < l; j++){
				const Poly32& polyModpj = vec_polys[j];
				int32_t ai_mod_pj = vec_polys[j][i];
				ai += (ai_mod_pj * Phat_invP[j]);
			}
			ai %= P;
			if (2*ai > P){
				ai -= P;
			}
			ai %= mod;
			if (2*ai > mod){
				ai -= mod;
			}
			ans[i] = (T) (ai.to<long>());
//			cout << "a"<<i<<" = " << ai << endl;
//
//
//			cout << "a"<<i<<".to<long> = " << ans[i] << endl;

		}
		return ans;
	}

	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	PolyCRT_FFT forward(const vector<T>& poly){
		PolyCRT_FFT out(l);
		vector<Poly32> polys = CRT_FFT::reduce(poly);

		for(int i = 0; i < l; i++){
			fft->to_fft(out[i], polys[i]);
		}

		return out;
	}

	template <typename T> // T is supposed to be Poly32 or Poly64
	vector<PolyCRT_FFT> forward_vector(const vector<T>& poly){
		int d = poly.size();
		vector<PolyCRT_FFT> out(d);

		for(int i = 0; i < d; i++){
			out[i] = forward(poly[i]);
		}

		return out;
	}

	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> backward(const PolyCRT_FFT& poly){
		
		vector<Poly32> inv_ffts(l);
		for(int i = 0; i < l; i++){
			fft->from_fft_core(poly[i]); // use fft core to avoid copying vector two times
			inv_ffts[i] = Poly32(N);
			for(int j = 0; j < N; j++){ // now copy vector already reducing it mod pi
				int64_t tmp =  (int64_t)(round(fft->in_array[j])); // this avoids overflows 
				inv_ffts[i][j] = (int32_t) (tmp % p[i]);
			}
		}
		vector<T> out = CRT_FFT::CRT<T>(inv_ffts);

		return out;
	}

	// Perform multiplication modulo X^N+1
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> multiply(const vector<T>& a, const vector<T>& b){
		PolyCRT_FFT va = forward(a);
		PolyCRT_FFT vb = forward(b);

		va *= vb;

		return backward<T>(va);
	}

	// compute the dot product of the subvectors a[begin ... begin+nelmnts-1]
	// and b[begin ... begin+nelmnts-1] modulo X^N+1
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product(const vector<vector<T> >& a, 
							  const vector<vector<T> >& b,
							  int begin,
							  int nelmnts){
		assert(a.size() == b.size());
		int d = a.size();
		PolyCRT_FFT va = forward(a[begin]);
		PolyCRT_FFT vb = forward(b[begin]);
		va *= vb;
		PolyCRT_FFT res(va);

		for(int i = begin+1; i < begin+nelmnts; i++){
			va = forward(a[i]);
			vb = forward(b[i]);
			// multiply entrywise
			va *= vb;
			// add all the entries
			res += va;
		}
		// now, res = va[begin]*vb[begin] + ... + va[begin+d-1]*vb[begin+d-1]

		return backward<T>(res);
	}




	// Compute the dot product a*b modulo X^N+1
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product(const vector<vector<T> >& a, const vector<vector<T> >& b){
		assert(a.size() == b.size());
		int d = a.size();

		// performs the inner product by chunks of d_max elements (or less)

		int rest = d % d_max;
		int n_partitions = (d - rest) / d_max;

		vector<T> res(N, 0);

		for(int i = 0; i < n_partitions; i++){
			res += inner_product(a, b, i * d_max, d_max);
		}
		if (0 != rest)
			res += inner_product(a, b, d - rest, rest);

		return res;
	}

	// compute the dot product of the subvectors a[begin ... begin+nelmnts-1]
	// and b[begin ... begin+nelmnts-1] modulo X^N+1
	// but assuming that the CRT-FFT of b has already been precomputed.
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product(const vector<vector<T> >& a, 
							  const vector<PolyCRT_FFT>& vb,
							  int begin,
							  int nelmnts){
		assert(a.size() == vb.size());
		int d = a.size();
		PolyCRT_FFT va = forward(a[begin]);
		va *= vb[begin];
		PolyCRT_FFT res(va);

		for(int i = begin+1; i < begin+nelmnts; i++){
			va = forward(a[i]);
			// multiply entrywise
			va *= vb[i];
			// add all the entries
			res += va;
		}
		// now, res = va[begin]*vb[begin] + ... + va[begin+d-1]*vb[begin+d-1]

		return backward<T>(res);
	}



	// 	Compute the dot product a*b modulo X^N+1, but receiving the CRT-FFT 
	// of b instead of b itself.
	//	This is useful when the same vector has to be multiplied by several 
	// different vectors.
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product(const vector<vector<T> >& a, const vector<PolyCRT_FFT>& vb)
	{
		assert(a.size() == vb.size());
		int d = a.size();

		// performs the inner product by chunks of d_max elements (or less)

		int rest = d % d_max;
		int n_partitions = (d - rest) / d_max;

		vector<T> res(N, 0);

		for(int i = 0; i < n_partitions; i++){
			res += inner_product(a, vb, i * d_max, d_max);
		}
		if (0 != rest)
			res += inner_product(a, vb, d - rest, rest);

		return res;
	}


	// compute the dot product of the subvectors a[begin ... begin+nelmnts-1]
	// and b[begin ... begin+nelmnts-1] modulo X^N+1
	// but assuming that the CRT-FFT of both a and b have already been precomputed.
	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product_full_precomp(const vector<PolyCRT_FFT >& va, 
							  const vector<PolyCRT_FFT>& vb,
							  int begin,
							  int nelmnts){
		assert(va.size() == vb.size());
		int d = va.size();
		PolyCRT_FFT res(va[begin]);
		res *= vb[begin];
		
		PolyCRT_FFT tmp(va[begin]);

		for(int i = begin+1; i < begin+nelmnts; i++){
			tmp = va[i];
			// multiply entrywise
			tmp *= vb[i];
			// add all the entries
			res += tmp;
		}
		// now, res = va[begin]*vb[begin] + ... + va[begin+d-1]*vb[begin+d-1]

		return backward<T>(res);
	}


	template <typename T> // T is supposed to be int32_t, int64_t, ZZ, etc
	vector<T> inner_product_full_precomp(const vector<PolyCRT_FFT>& va, const vector<PolyCRT_FFT>& vb)
	{
		assert(va.size() == vb.size());
		int d = va.size();

		// performs the inner product by chunks of d_max elements (or less)

		int rest = d % d_max;
		int n_partitions = (d - rest) / d_max;

		vector<T> res(N, 0);

		for(int i = 0; i < n_partitions; i++){
			res += inner_product_full_precomp<T>(va, vb, i * d_max, d_max);
		}
		if (0 != rest)
			res += inner_product_full_precomp<T>(va, vb, d - rest, rest);

		return res;
	}


};

std::ostream& operator<< (std::ostream &out, const CRT_FFT& u);


#endif

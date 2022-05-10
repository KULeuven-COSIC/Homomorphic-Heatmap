
#ifndef __PolynomialMultiplier__
#define __PolynomialMultiplier__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

#include "NTT.h"
#include "utils.h"

#include <vector>

#define PrecomputedVecNTT std::vector< NTL::vec_ZZ >


/**
 * 		Let R := Z[x] / <x^N+1>, with N being a power of two.
 *		This class is used to perform efficient inner products of vectors 
 *	over R.
 *
 *		It uses the class NTT two perform efficient negacyclic convolutions,
 *	which are equivalent to multiplication of polynomials of R.
 **/
class PolynomialMultiplier {
	

	public:

		NTT* ntt;

		long int d;  // dimension of vectors that will be multiplied

		NTL::ZZ psi; // primitive 2N-th root of unity on Z/pZ (order(psi) = 2N in Z/pZ)
		NTL::ZZ inv_psi; // psi^-1 in Z/pZ
		NTL::vec_ZZ powers_psi;  // precomputed psi^i % p for 0 <= i < N
		NTL::vec_ZZ powers_inv_psi; // precomputed psi^-i % p for 0 <= i < N
		
		// auxiliary vectors used to store a polynomial with the result and one
		// of the operands during the computation of the inner products.
		NTL::vec_ZZ vres; 
		NTL::vec_ZZ vai;


	/** 
	 *		logN is the logarithm of N in base 2
	 *		d is the dimension of the vectors that will be multiplied
	 *		b0 is a bound to the log infinity norm of one of the vectors
	 *		b1 is a bound to the log infinity norm of the other vector
	 */
	PolynomialMultiplier(int logN, int d, int b0, int b1);

	~PolynomialMultiplier();

	/**
	 * 		Let N be the dimension of the NTT, sigma the bit reverse order
	 * 	permutaion precomputed in NTT, and o be the entrywise vector product modulo p.
	 * 		This function copies a to an N-dimensional vector u, and returns
	 * 	sigma(u o powers).
	 */
	NTL::vec_ZZ apply_bit_reverse_and_powers_to_poly(
						const NTL::ZZX& a, 
						const NTL::vec_ZZ& powers);


	/** Precompute the weighted NTT of a (multiplied by powers of psi) to speedup
	 * future products by a.
	 */
	NTL::vec_ZZ precompute_NTT(const NTL::ZZX& a);

	/** Precompute the weighted NTT of all entries of v to speedup future 
	 * products by v.
	 */
	PrecomputedVecNTT precompute_NTT(const std::vector<NTL::ZZX >& v);


	/**
	 *		Receives a polynomial a and a precomputed NTT of a polynomial b
	 * and returns a polynomial equal to a*b mod x^N+1.	 					*/
	NTL::ZZX multiply(const NTL::ZZX& a, const NTL::vec_ZZ& vb);

	/**
	 *		Receives a vector of polynomials va and 
	 *	vector of precomputed NTTs and returns a polynomial equal
	 *	to inner product modulo x^N+1.	 					*/
	NTL::ZZX multiply(const std::vector<NTL::ZZX>& va, const PrecomputedVecNTT& vb);

};

#endif

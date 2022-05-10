#ifndef __SIMPLE_FV_SCHEME__
#define __SIMPLE_FV_SCHEME___

#include <vector>
#include "flint/fmpz_polyxx.h"
#include <random>
#include <chrono>

#include "utils.h"

using namespace flint;

class Ciphertext {
	public:
		fmpz_polyxx a;
		fmpz_polyxx b;


	Ciphertext();

	Ciphertext(const fmpz_polyxx& _a, const fmpz_polyxx& _b);

};

class NonRelinCiphertext {
	public:
		fmpz_polyxx a;
		fmpz_polyxx b;
		fmpz_polyxx c;


	NonRelinCiphertext();

	NonRelinCiphertext(const fmpz_polyxx& _a, const fmpz_polyxx& _b, const fmpz_polyxx& _c);
};


//	Plaintext space is the ring Z_t[x]/<x^N + 1>.
class FV {

	public:
		long logq;
		fmpzxx q;
		long N; // The polynomial ring is R = ZZ[x] / <x^N + 1>
		long logN;
		fmpzxx B; // base in which we decompose ciphertexts during multiplication
		long logB; // logarithm of B in base 2
		long l; // logarightm of q in base B
		double sigma; 

		bool bin_sk;

		fmpz_polyxx fmod; // fmod = x^N + 1, defines the degree of the plaintext space
		fmpzxx t;  // defines the maximum value of the coefficients of plaintexts
		fmpzxx Delta; // q/t

		long alpha; // value used to the automorphism in the heatmap application

    	std::default_random_engine rand_gen;
		std::normal_distribution<double> gaus_dist; 

	//private:
		fmpz_polyxx s; // degree-(N-1) random polynomial with binary coefficients
		

		std::vector<fmpz_polyxx> words_vec; // vector to decompose polynomials during key-switching
		std::vector<Ciphertext> rlk; // relinearization key

		std::vector<Ciphertext> autk; // key-switching key from sk(X^alpha) to sk


	FV();

	FV(long logN, long logq, double sigma, long logB, long t, long alpha, bool bin_sk=true);

	void gen_sk();
	
	void gen_relin_key();
	
	// Generates key-switching key to allow the automorphism X --> X^alpha
	void gen_automorphism_key();

	Ciphertext enc(const fmpz_polyxx& m);
	
	Ciphertext enc_on_exponent(long m);

	fmpz_polyxx dec(const Ciphertext& c);

	fmpz_polyxx dec(const NonRelinCiphertext& c);


	NonRelinCiphertext mult_without_relin(const Ciphertext& c0, const Ciphertext& c1);
	// Multiply and relinearize	
	Ciphertext mult(const Ciphertext& c0, const Ciphertext& c1);

	double get_noise(const Ciphertext& c, const fmpz_polyxx& msg);
	
	double get_noise(const NonRelinCiphertext& c, const fmpz_polyxx& msg);

	// Apply automorphism X --> X^alpha to c.
	// If key_swt is true, then apply key-switching from sk(X^alpha) back to sk
	Ciphertext apply_automorphism(const Ciphertext& c, bool key_swt=true);


	// auxiliary functions
	fmpz_polyxx sample_noise();
	
};
	
std::ostream& operator<< (std::ostream &out, const FV& he);

#endif


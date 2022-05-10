#ifndef __HOMOMORPHIC_HEATMAP__
#define __HOMOMORPHIC_HEATMAP__

#include <vector>

#include "flint/fmpz_polyxx.h"

#include "utils.h"
#include "FV.h"
#include "CRT_FFT.h"

using namespace flint;
using namespace std;

typedef std::vector<Ciphertext> vector_ciphertext;

// This class model the ciphertext we obtain after multiplying an RLWE ciphertext
// by a test vector defined with powers of x.
class MixedCiphertext {
	public:
	std::vector<fmpz_polyxx> a;
	fmpz_polyxx b; // b = a * vec(sk) + e + Delta * m

	MixedCiphertext(const std::vector<fmpz_polyxx>& _a, const fmpz_polyxx& _b);
};


//	Plaintext space is the ring Z_t[x]/<x^N + 1>.
class HomomorphicHeatmap {

	public:
		// Assuming a heatmap with x coordinates from 0 to x_max - 1, y coordinates
		// from 0 to y_max -1, and cells with dimensions base_cell x height_cell.
		// Thus, the function that maps each point (x, y) to its cell is
		// 			f(x) * alpha + g(y)
		// 	with 
		// 		f(x) = floor(x / base_cell)
		// 		g(y) = floor(y / height_cell)
		// 		alpha = ceil(y_max / base_cell)
		long x_max;
		long y_max;
		long base_cell;
		long height_cell;



		long alpha; 

		// Parameters for FV scheme that will compute X^f(x)
		long log_qx;
		fmpzxx qx;
		long log_Px;
		fmpzxx Px;
		long Nx; // The polynomial ring is R = ZZ[x] / <x^N + 1>
				// We need Nx >= floor((x_max - 1) / base_cell)
		int log_Nx;
		long k_x; // number of sets in the partition of the domain [[0, x_max[[
		double sigma_x; 

		// Parameters for FV scheme that will compute X^g(y)
		long log_qy;
		fmpzxx qy;
		long log_Py;
		fmpzxx Py;
		long Ny; // The polynomial ring is R = ZZ[x] / <x^N + 1>
				// We need Ny >= floor((y_max - 1) / height_cell)
		int log_Ny;
		long k_y; // number of sets in the partition of the domain [[0, x_max[[
		double sigma_y; 


		std::vector<long> func_f; // k_x * Nx dimensional vector
		std::vector<long> func_g; // k_y * Ny dimensional vector

		long img_size_f;
		long img_size_g;
		long final_img_size;

		int Nbar; // The polynomial ring is R = ZZ[x] / <x^N + 1>
				  // We need Nbar >= ceil(x_max / base_cell) * alpha + ceil(y_max / height_cell)
				  //               > max(f(x)) * alpha + max(g(y)) 
		long log_Nbar;


		fmpzxx t;  // number of elements that will be counted
		
		fmpzxx Delta;  // q / t rounded

		FV fvx;
		FV fvy;
		FV fvbar;

		CRT_FFT* crt_fft; // to compute FFTs, multiply on R, etc

		std::vector<PolyCRT_FFT> ffk_a_x; // format-fixing key. Nx-dimensional vector over Rbar
		std::vector<PolyCRT_FFT> ffk_b_x; // format-fixing key. Nx-dimensional vector over Rbar
		std::vector<PolyCRT_FFT> ffk_a_y; // format-fixing key. Ny-dimensional vector over Rbar
		std::vector<PolyCRT_FFT> ffk_b_y; // format-fixing key. Ny-dimensional vector over Rbar


		vector<vector<fmpzxx> > res_vec_x; // Auxiliar ceil(x_max / base_cell) x Nx matrix
		vector<vector<fmpzxx> > res_vec_y; // Auxiliar ceil(y_max / height_cell) x Ny matrix

	//private:

	HomomorphicHeatmap();

	HomomorphicHeatmap(long x_max, long y_max, long base, long height, long log_Nx, long log_Ny, long log_qx, long log_qy, long t);

//	~HomomorphicHeatmap();

	// Returns a vector ciphertext composed of k ciphertexts 
	// where one of them encrypts a power of X and all
	// the others encrypt zero.
	// char x_or_y determines if fvx or fvy are used to encrypt.
	vector_ciphertext enc_high_level(long m, char x_or_y);

	// Returns a vector of "vector ciphertexts", each one composed by k standard RLWE
	// ciphertexts. Each vector ciphertext encrypts one entry of ms 
	std::vector< vector_ciphertext > enc_high_level(const std::vector<long>& ms, char x_or_y);


	/* 	Generates the two format-fixing keys:
	 *  	- the first one encrypts coef_vec(fvy.sk) under fvbar.sk
	 *  	 modulus P*q and X^Nbar+1, i.e., we have RLWE samples encrypting an
	 *  	 integer vector;
	 *  	- the second one encrypts coef_vec(fvx.sk) under fvbar.sk(X^beta)
	 *  	where beta * alpha = 1 mod 2*N.
	 *
	 *  This beta is needed to skip the key-switching key that we would need
	 * after applying the autormorphism to get enc(X^(alpha * f(x)). 
	 *
	 *  With this, the size of the key is 2*N * bar{N} * log(P*q). */
	void gen_format_fixing_keys();


	/* Receives an N-dimensional integer vector u whose entries are in [[0, img_size-1]]
	* and a polynomial b of degree <= N. 
	* Let w be the vector 
	* 		(X^u[0], ..., X^u[N-1]),
	* i.e., powers of X defined by u,
	* and let v be the vector of coefficients of b.
	* This function returns the polynomial w*v */
	fmpz_polyxx mult_powers_X_by_vec(const std::vector<long>& u, 
								  const fmpz_polyxx& b,
								  long N,
								  long img_size);


	/* Receives an N-dimensional integer vector u whose entries are in [[0, img_size-1]]
	* and a polynomial a of degree smaller than N.
	*   Let w be the vector (x^u[0], ..., x^u[d-1]), so powers of x defined by u.
	*   Let A be the NxN circulant matrix of a modulo x^d + 1.
	*   Set res equal to the vector of polynomials w*A */
	void mult_powers_X_by_circ_mat_by_diag(
												vector<fmpz_polyxx>& res,
												const std::vector<long>& u, 
											  	const fmpz_polyxx& a,
												long img_size,
												const fmpzxx& q,
												char x_or_y
											  );

	/* Receives a RLWE ciphertext (a, b) defined mod (x^N + 1, q) and an 
	 * N-dimensional vector representing one partition of the function */
	MixedCiphertext mult_by_test_vector(const Ciphertext& c,
										const std::vector<long>& func,
										long N,
										long img_size,
										const fmpzxx& q,
										char x_or_y
										);


//	std::vector<long> embed_func(const vector<long>& func);


	/* cs: k ciphertexts as generated by enc_high_level
	 * char x_or_y: if 'x', then apply f, if 'y', then apply g
	**/
	MixedCiphertext apply_func_to_exponent(const vector_ciphertext & cs,
									  	  char x_or_y);

	/* Receives an "anomalous" RLWE ciphertext c in R^(Nx+1) or R^(Ny + 1), 
	 * as returned by apply_func_to_exponent.
	* 	Fix the format of c transforming it in a normal RLWE ciphertext mod x^Nbar+1
	* using the key ffk_x or ffk_y depending if the char ffk is 'x' or 'y'.
	**/
	Ciphertext apply_format_fixing_key(const MixedCiphertext& c, char ffk);


	/* cs_x: k ciphertexts as generated by enc_high_level.
	 * cs_y: k ciphertexts as generated by enc_high_level.
	 * Assuming that (cs_x, cs_y) encrypts a point (x, y), this function
	 * returns an RLWE encryption of X^(f(x) * alpha + g(y))
	 */
	NonRelinCiphertext compute_cell(
							const vector_ciphertext & cs_x,
							const vector_ciphertext & cs_y
							);
	

	/* cs_x: vector of vector ciphertexts encrypting the x coordinates.
	 * cs_y: vector of vector ciphertexts encrypting the y coordinates.
	 * Both cs_x and cs_y must be obtained from the function enc_high_level.
	 *
	 * Returns an RLWE ciphertext encrypting a polynomial that represents the 
	 * heatmap (each coefficient a_i of the polynomial represents how many points
	 * lie in the i-th cell).
	 */
	NonRelinCiphertext compute_heatmap(
								const std::vector< vector_ciphertext >& cs_x,
								const std::vector< vector_ciphertext >& cs_y
								);
};

std::ostream& operator<< (std::ostream &out, const HomomorphicHeatmap& hhm);

#endif


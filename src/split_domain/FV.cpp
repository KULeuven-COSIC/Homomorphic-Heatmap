#include "FV.h"

using namespace flint;

Ciphertext::Ciphertext() {  };

Ciphertext::Ciphertext(const fmpz_polyxx& _a, const fmpz_polyxx& _b) : a(_a), b(_b) {  };

NonRelinCiphertext::NonRelinCiphertext() {  };

NonRelinCiphertext::NonRelinCiphertext(const fmpz_polyxx& _a, 
										const fmpz_polyxx& _b, 
										const fmpz_polyxx& _c)
										: a(_a), b(_b), c(_c) {  };

FV::FV() {

}

FV::FV(long logN, long logq, double sigma, long logB, long t, long automorphism_exp, bool bin_sk) 
	: gaus_dist(0, sigma / 2.5066) // sqrt(2*pi)
{
	this->logN = logN;
	this->logq = logq;
	this->logB = logB;
	this->B = fmpzxx(1) << logB;
	this->sigma = sigma;
	this->bin_sk = bin_sk;

    unsigned seed = chrono::steady_clock::now().time_since_epoch().count(); 
    this->rand_gen = std::default_random_engine(seed);

	if (0 != automorphism_exp && 0 == automorphism_exp % 2){
		cout << "ERROR: automorphism_exp must be zero, which means "
			 << "that no automorphism functionality is desired, or an "
			 << "odd number." << endl;
		exit(15);
	}
	this->alpha = automorphism_exp;

	this->N = 1 << logN;
	// fmod = x^N + 1
	fmod.set_coeff(N, 1); 
	fmod.set_coeff(0, 1);

	this->q = fmpzxx(1);
	this->q = q << logq;

	this->t = t;
	this->Delta = rounded_division(q, fmpzxx(t));

	this->l = ceil(logq / (double)logB) + 1; // logarightm of q in base B
	this->words_vec = vector<fmpz_polyxx>(l);

	this->gen_sk();

	this->gen_relin_key();

	this->gen_automorphism_key();
}

void FV::gen_sk(){
	if(bin_sk){
		this->s = random_poly(N, 2);
		reduce_poly_mod(s, fmpzxx(2));
	}else{
		this->s = random_poly(N, logq);
		reduce_poly_mod(s, q);
	}
}
	
void FV::gen_relin_key() {
	fmpz_polyxx pow_s_square(s * s);
	reduce_2N_to_N(pow_s_square, N);
	for(long int i = 0; i < l; i++){
		Ciphertext ci = enc(fmpz_polyxx(0)); // b = -a*s + e mod q
		ci.b += pow_s_square;
    	reduce_poly_mod(ci.b, q); // b = -a*s - e + B^i * s^2 mod q
		rlk.push_back(ci);
		pow_s_square *= B;
	}
}
	
void FV::gen_automorphism_key() {
	fmpz_polyxx pow_s_alpha = eval_poly(s, alpha, N);
	for(long int i = 0; i < l; i++){
		Ciphertext ci = enc(fmpz_polyxx(0)); // b = -a*s + e mod q
		ci.b += pow_s_alpha;
    	reduce_poly_mod(ci.b, q); // b = -a*s - e + B^i * s(X^alpha) mod q
		autk.push_back(ci);
		pow_s_alpha *= B;
	}
}

Ciphertext FV::enc(const fmpz_polyxx& m) {
	fmpz_polyxx a = random_poly(N, logq);
    fmpz_polyxx e = this->sample_noise();
    fmpz_polyxx b(a*s + e + Delta * m);
	b %= fmod;
    reduce_poly_mod(b, q);
    return Ciphertext(a, b);
}

Ciphertext FV::enc_on_exponent(long m){
	fmpz_polyxx msg; 
	msg.set_coeff(m, 1); // msg = x^m
	return enc(msg);
}


fmpz_polyxx FV::dec(const Ciphertext& c) {
    fmpz_polyxx noisy_msg(c.b - c.a * s);
	noisy_msg %= fmod;
    reduce_poly_mod(noisy_msg, q);
	noisy_msg *= t;
    fmpz_polyxx m = rounded_division(noisy_msg, q);
    reduce_poly_mod(m, t, false);
	return m;
}


fmpz_polyxx FV::dec(const NonRelinCiphertext& c) {
	fmpz_polyxx noisy_msg(c.c - c.b * s + c.a * s * s);
	noisy_msg %= fmod;
    reduce_poly_mod(noisy_msg, q);
	noisy_msg *= t;
    fmpz_polyxx m = rounded_division(noisy_msg, q);
    reduce_poly_mod(m, t, false);
	return m;
}

NonRelinCiphertext FV::mult_without_relin(const Ciphertext& c0, const Ciphertext& c1) {
	fmpz_polyxx a0( t*(c0.a * c1.a) );
	fmpz_polyxx a1( t*(c0.a * c1.b + c0.b * c1.a) );
	fmpz_polyxx a2( t * c0.b * c1.b );
	
	reduce_2N_to_N(a0, N);
	reduce_2N_to_N(a1, N);
	reduce_2N_to_N(a2, N);

	a0 = rounded_division(a0, q);
	a1 = rounded_division(a1, q);
	a2 = rounded_division(a2, q);

	reduce_poly_mod(a0, q);
	reduce_poly_mod(a1, q);
	reduce_poly_mod(a2, q);

	return NonRelinCiphertext(a0, a1, a2);
}

Ciphertext FV::mult(const Ciphertext& c0, const Ciphertext& c1) { 

	NonRelinCiphertext c = mult_without_relin(c0, c1);

	fmpz_polyxx& a0 = c.a;
	fmpz_polyxx& a1 = c.b;
	fmpz_polyxx& a2 = c.c;

	// now, relinearize:
	decompose_ZZX(words_vec, a0, N, l, B);
	for(long i = 0; i < l; i++){
		a2 += (rlk[i].b * words_vec[i]);
		reduce_2N_to_N(a2, N);
		a1 += (rlk[i].a * words_vec[i]);
		reduce_2N_to_N(a1, N);
	}
	reduce_poly_mod(a1, q);
	reduce_poly_mod(a2, q);

	return Ciphertext(a1, a2);
}

Ciphertext FV::apply_automorphism(const Ciphertext& c, bool key_swt){
	fmpz_polyxx a = eval_poly(c.a, alpha, N);
	fmpz_polyxx b = eval_poly(c.b, alpha, N);

	if (key_swt){
		// now, key-switch:
		decompose_ZZX(words_vec, a, N, l, B);
		a = fmpz_polyxx(0);
		for(long i = 0; i < l; i++){
			b -= autk[i].b * words_vec[i];
			reduce_2N_to_N(b, N);
			a -= autk[i].a * words_vec[i];
			reduce_2N_to_N(a, N);
		}
		reduce_poly_mod(a, q);
		reduce_poly_mod(b, q);
	}

	return Ciphertext(a, b);
}

double FV::get_noise(const Ciphertext& c, const fmpz_polyxx& msg) {
    fmpz_polyxx noise(c.b - c.a * s - Delta * msg);
	noise %= fmod;
    reduce_poly_mod(noise, q);
	fmpzxx i_n = max_norm(noise);
	if (i_n.is_zero())
		return -1.0;
	return flint::dlog(i_n) / log(2);
}
	
double FV::get_noise(const NonRelinCiphertext& c, const fmpz_polyxx& msg){
	fmpz_polyxx noise(c.c - c.b * s + c.a * s * s - Delta*msg);
	noise %= fmod;
    reduce_poly_mod(noise, q);
	fmpzxx i_n = max_norm(noise);
	if (i_n.is_zero())
		return -1.0;
	return flint::dlog(i_n) / log(2);
}

fmpz_polyxx FV::sample_noise(){ 
	fmpz_polyxx e;
	for (long i = N-1; i >= 0; i--){
		long int r = (long int) this->gaus_dist(rand_gen);
		fmpzxx ei(r);
		e.set_coeff(i, ei);
	}
	return e;
}
	
std::ostream& operator<< (std::ostream &out, const FV& he) {  
	out << "FV: {" 
	   << "N: 2^" << he.logN
	   << ", q: 2^" << he.logq
	   << ", sigma: " << he.sigma
	   << ", t: " << he.t
	   << ", B: " << he.B
	   << ", l: " << he.l
	   << ", alpha: " << he.alpha
	   << ", " << (he.bin_sk? "bin sk" : "uniform sk")
	   << "}";
	return out;
}


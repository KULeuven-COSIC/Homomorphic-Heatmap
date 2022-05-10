#include "CRT_FFT.h"
#include "utils.h"
#include <cassert>

using namespace flint;

int primes15bits[] = {32771,32779,32783,32789,32797,32801,32803,32831,32833,
					 32839,32843,32869,32887,32909,32911,32917,32933,32939, 
					 32941,32957,32969,32971,32983,32987,32993,32999,33013,
					 33023,33029,33037,33049,33053,33071,33073,33083,33091,
					 33107,33113,33119,33149,33151,33161,33179,33181,33191,
					 33199,33203,33211,33223,33247}; // first 50 primes with 15-bits



int primes20bits[] = {1048583, 1048589, 1048601, 1048609, 1048613, 
					  1048627, 1048633, 1048661, 1048681, 1048703, 
					  1048709, 1048717, 1048721, 1048759, 1048783, 
					  1048793, 1048799, 1048807, 1048829, 1048837, 
					  1048847, 1048867, 1048877, 1048889, 1048891, 
					  1048897, 1048909, 1048919, 1048963, 1048991}; // first 30 primes with 20-bits


CRT_FFT::CRT_FFT(int logN, int bl, int br, int d){
	N = 1<<logN;
	this->logN = logN;
	this->bl = bl;
	this->bl = br;
	this->d_max = 80; // XXX: maybe find a way of computing this automatically
					 // increasing this value increases the changes of getting an approximation error
					 // decreasing it slows down the inner products
	
	this->fft = new FFT_engine(N);

	int needed_bits = logN / 2 + bl + br + log(sqrt(d))/log(2);

	this->l = ceil(needed_bits / 20.0); // assuming each prime has 20 bits

	this->p = vector<int>(l);
	this->invP = vector<int>(l);

	P = fmpzxx(1);

	for(int i = 0; i < l; i++){
		p[i] = ::primes20bits[i];
		P *= p[i];
	}

	Phat = vector<flint::fmpzxx>(l); // Phat[i] = P / pi
	for(int i = 0; i < l; i++){
		Phat[i] = P / p[i];
	}

	compute_inverses_crt(invP, p);


	Phat_invP = vector<flint::fmpzxx>(l); // Phat[i] = (P / pi) * ((P/pi)^-1 mod pi) mod P
	for(int i = 0; i < l; i++){
		Phat_invP[i] = (Phat[i] * invP[i]) % P;
	}
}

CRT_FFT::~CRT_FFT(){
	delete this->fft;
}

std::ostream& operator<< (std::ostream &out, const CRT_FFT& u){
	out << "CRT_FFT: {" 
	   << "N: 2^" << u.logN
	   << ", l: " << u.l
	   << ", primes: " << u.p
	   << "}";
	return out;

}




#include "PolynomialMultiplier.h"

using namespace NTL;
using namespace std;

PolynomialMultiplier::PolynomialMultiplier(int logN, int d, int b0, int b1){

	this->d = d;

	long int N = 1 << logN;
	NTL::ZZ p = find_modulus(b0 + ceil(log(d)/log(2)), b1, 2*N);
	this->psi = find_Nth_primitive_root(2*N, p);
	this->inv_psi = InvMod(psi, p);
	NTL::ZZ omega = psi*psi % p;

	this->ntt = new NTT(logN, p, omega);

	this->powers_psi.SetLength(N);
	this->powers_inv_psi.SetLength(N);
	powers_psi[0] = ZZ(1);
	powers_inv_psi[0] = ZZ(1);

	// precompute powers of psi and psi^-1
	for (long int i = 1; i < N; i++){
		powers_psi[i] = (powers_psi[i-1] * psi) % p;
		powers_inv_psi[i] = (powers_inv_psi[i-1] * inv_psi) % p;
	}

	this->vres.SetLength(N);
	this->vai.SetLength(N);
	for (long int i = 0; i < N; i++){
		vres[i] = ZZ(0);
		vai[i] = ZZ(0);
	}

}
	
PolynomialMultiplier::~PolynomialMultiplier(){
	free(this->ntt);
}

vec_ZZ PolynomialMultiplier::apply_bit_reverse_and_powers_to_poly(
					const ZZX& a, 
					const vec_ZZ& powers){

	ZZ& p = ntt->p;
	long int N = ntt->N;
	vector<long int>& bit_rev_perm = ntt->bit_rev_perm; 
	vec_ZZ va; va.SetLength(N);
	for (long int i = 0; i < N; i++) {
		long int j = bit_rev_perm[i];
		if (i == j)
			va[i] = coeff(a, i) * powers[i] % p;
		else if (j > i){
			va[i] = coeff(a, j) * powers[j] % p;
			va[j] = coeff(a, i) * powers[i] % p;
		}
	}
	return va;
}


// NTT( (psi^0, psi^1, ..., psi^(N-1)) x (a0, a1, ..., a_(N-1)) )
NTL::vec_ZZ PolynomialMultiplier::precompute_NTT(const NTL::ZZX& a){
	vec_ZZ va = apply_bit_reverse_and_powers_to_poly(a, powers_psi);
	ntt->low_level_transform(va, ntt->powers_omega);
	return va;
}

PrecomputedVecNTT PolynomialMultiplier::precompute_NTT(const std::vector<NTL::ZZX >& v){
	vector< vec_ZZ > va = vector< vec_ZZ >(v.size());
	for (long int i = 0; i < v.size(); i++)
		va[i] = precompute_NTT(v[i]);
	return va;
}


NTL::ZZX PolynomialMultiplier::multiply(const NTL::ZZX& a, const NTL::vec_ZZ& vb){

	ZZ& p = ntt->p;
	// NTT(powers_psi o a) where o is entrywise multiplication
	vres = apply_bit_reverse_and_powers_to_poly(a, powers_psi); 
	ntt->low_level_transform(vres, ntt->powers_omega);

	// multiply
	coordinatewise_prod(vres, vb, vres, p);

	// NTT^-1(powers_psi^-1 o prod) mod p
	ntt->inv_transform(vres);
	coordinatewise_prod(vres, powers_inv_psi, vres);
	centralized_mod(vres, p);

	return vec_ZZ_to_poly(vres);
}
	
ZZX PolynomialMultiplier::multiply(const vector<ZZX>& va, const PrecomputedVecNTT& vb){
	ZZ& p = ntt->p;
	vres = apply_bit_reverse_and_powers_to_poly(va[0], powers_psi);
	ntt->low_level_transform(vres, ntt->powers_omega);
	coordinatewise_prod(vres, vb[0], vres, p);

	for (long int i = 1; i < d; i++){
		vai = apply_bit_reverse_and_powers_to_poly(va[i], powers_psi);
		ntt->low_level_transform(vai, ntt->powers_omega);
		
		coordinatewise_prod(vai, vb[i], vai, p);
		vres += vai;
	}

	ntt->inv_transform(vres);
	apply_powers(vres, inv_psi, p);
	centralized_mod(vres, p);

	return vec_ZZ_to_poly(vres);
}


#include "HomomorphicHeatmap.h"

using namespace std;

using namespace std::chrono; // to measure execution times

#define VERBOSE true

#if VERBOSE
// To measure run time
high_resolution_clock::time_point _t_before;
high_resolution_clock::time_point _t_after;
double _t_duration;
#endif

MixedCiphertext::MixedCiphertext(const std::vector<fmpz_polyxx>& _a, const fmpz_polyxx& _b) {
	this->a = vector<fmpz_polyxx>(_a.size());
	for(long i = 0; i < _a.size(); i++)
		a[i] = _a[i];
	this->b = _b;
}

HomomorphicHeatmap::HomomorphicHeatmap(long x_max, long y_max, long base, long height, 
										long log_Nx, long log_Ny, 
										long log_qx, long log_qy, 
										long t){
	this->x_max = x_max;
	this->y_max = y_max;
	this->base_cell = base;
	this->height_cell = height;

	this->alpha = ceil(y_max / height); // the final function is f(x) * alpha + g(y)
	if (0 == alpha % 2)
		alpha--; // alpha must be odd
	
	this->log_Nx = log_Nx;
	this->Nx = 1 << log_Nx;
	if (Nx <= floor((x_max - 1) / base_cell)){
		cout << "ERROR: Nx = " << Nx
			 << " is not big enough to cover the image of f, which is "
			 << "[[0, " << floor((x_max-1) / base_cell) << "]]" << endl;
		exit(1);
	}
	this->k_x = ceil(x_max / Nx); // domain is partitioned as P_1 U P_2 U ... U P_k
	if (0 == k_x)
		k_x = 1;


	this->log_Ny = log_Ny;
	this->Ny = 1 << log_Ny;
	if (Ny <= floor((y_max-1) / height_cell)){
		cout << "ERROR: Ny = " << Ny
			 << " is not big enough to cover the image of f, which is "
			 << "[[0, " << floor((y_max-1) / height_cell) << "]]" << endl;
		exit(2);
	}
	this->k_y = ceil(y_max / Ny); // domain is partitioned as P_1 U P_2 U ... U P_k
	if (0 == k_y)
		k_y = 1;
	
	this->log_qx = log_qx;
	this->log_Px = log_qx;
	this->qx = fmpzxx(1) << log_qx;
	this->Px = fmpzxx(1) << log_Px;

	this->log_qy = log_qy;
	this->log_Py = log_qy;
	this->qy = fmpzxx(1) << log_qy;
	this->Py = fmpzxx(1) << log_Py;

	this->t = t;


	this->func_f = vector<long>(k_x * Nx);
	this->func_g = vector<long>(k_y * Ny);
	for(long i = 0; i < x_max; i++)
		func_f[i] = floor(i / base_cell);
	for(long i = x_max; i < k_x * Nx; i++)
		func_f[i] = 0;

	for(long i = 0; i < y_max; i++)
		func_g[i] = floor(i / height_cell);
	for(long i = y_max; i < k_y * Ny; i++)
		func_g[i] = 0;


//	this->Delta = rounded_division(q, fmpzxx(t));

	this->final_img_size = floor((x_max-1) / base_cell) * alpha + floor((y_max-1) / height_cell);
	cout << "final_img_size = " << final_img_size << endl;
	cout << "log(final_img_size) = " << log(final_img_size) / log(2) << endl;
	this->log_Nbar = ceil(log(final_img_size+1) / log(2));
	if (log_Nbar < 12)
		this->log_Nbar = 12;
	this->Nbar = 1 << log_Nbar;

	// Initialize object of CRT-FFT to precompute FFTs of the fix-format keys
	crt_fft = new CRT_FFT(log_Nbar, 32, 64, 80); // there is no point in choosing last parameter larger than 80 because that is the
												 // d_max in the CRT_FFT class
	
	cout << *crt_fft << endl;

	double sigma_x = 3000;
	double sigma_y = 3000;
	double sigma_bar = 3.2;

	this->fvx = FV(log_Nx, log_qx, sigma_x, 4, t, /*alpha=*/0, /*bin_sk=*/false);
	cout << fvx << endl;
	
	this->fvy = FV(log_Ny, log_qy, sigma_y, 4, t, /*alpha=*/0,/*bin_sk=*/false);
	cout << fvy << endl;

	// XXX: the following line assumes log_qx = log_qy
	this->fvbar = FV(log_Nbar, log_qx, sigma_bar, 2, t, alpha);
	cout << fvbar << endl;

	this->img_size_f = floor((x_max - 1) / base_cell) + 1;
	this->res_vec_x = vector<vector<fmpzxx> >(img_size_f); 
	for(long i = 0; i < img_size_f; i++){
		res_vec_x[i] = vector<fmpzxx>(Nx);
	}

	this->img_size_g = floor((y_max - 1) / height_cell) + 1;
	this->res_vec_y = vector<vector<fmpzxx> >(img_size_g);
	for(long i = 0; i < img_size_g; i++){
		res_vec_y[i] = vector<fmpzxx>(Ny);
	}


	this->gen_format_fixing_keys(); // init ffk_x and ffk_y
}

std::vector<Ciphertext> HomomorphicHeatmap::enc_high_level(long m, char x_or_y){
	long k = ('x' == x_or_y ? k_x : k_y);
	long N = ('x' == x_or_y ? Nx : Ny);
	FV& fv = ('x' == x_or_y ? fvx : fvy);

    vector<Ciphertext> cs(k);
    long j = m % N;
    long i = (m - j) / N; // m belongs to the i-th partition of the domain
	fmpz_polyxx zero(0);

	for(long v = 0; v < i; v++)
		cs[v] = fv.enc(zero);

    cs[i] = fv.enc_on_exponent(j);

	for(long v = i+1; v < k; v++)
		cs[v] = fv.enc(zero);

    return cs;
}

std::vector< vector_ciphertext > HomomorphicHeatmap::enc_high_level(const std::vector<long>& ms, char x_or_y){
	long n = ms.size();
	std::vector< vector_ciphertext > ctxts(n);
	for(long i = 0; i < n; i++)
		ctxts[i] = enc_high_level(ms[i], x_or_y);
	return ctxts;
}

// Sets ffk as a format-fixing key from fv.sk to fvbar.sk(X^beta)
void gen_format_fixing_key(
							vector<PolyCRT_FFT>& ffk_a,
							vector<PolyCRT_FFT>& ffk_b,
							FV& fv, 
							FV& fvbar,
							long logP, 
							long logq, 
							long beta,
							CRT_FFT* crt_fft)
{

	fmpzxx P(fmpzxx(1) << logP);
	fmpzxx q(fmpzxx(1) << logq);
	fmpz_polyxx Ps(P * fv.s);
	fmpz_polyxx sbar(fvbar.s);

	if (1 != beta)
		sbar = eval_poly(sbar, beta, fvbar.N); // sbar = fv.sk(X^beta)

	fmpz_polyxx a, b, e;
	int logPq = logP + logq;
	fmpzxx Pq(P * q);

	ffk_a = vector<PolyCRT_FFT>(fv.N); // alocate space for N ciphertexts
	ffk_b = vector<PolyCRT_FFT>(fv.N); // alocate space for N ciphertexts
	for(long i = 0; i < fv.N; i++){

		a = random_poly(fvbar.N, logPq);
		reduce_poly_mod(a, Pq);
		e = fvbar.sample_noise();

		b = a * sbar;
		b += e;
		b.coeff(0) += Ps.get_coeff(i);
		reduce_2N_to_N(b, fvbar.N);  // b %= fvbar.fmod;
		reduce_poly_mod(b, Pq);
		
		a *= -1; // when ffk is applied, we have to compute -c.a * ffk_a, so multiply by -1 here already
		ffk_a[i] = crt_fft->forward(polyFlint_to_Poly64(a));
		ffk_b[i] = crt_fft->forward(polyFlint_to_Poly64(b));
		if (0 == i % 100)
			cout << "i = " << i << endl;
	}
	cout << "finished gen_format_fixing_key" << endl;
}


void HomomorphicHeatmap::gen_format_fixing_keys(){
	long twoN = 2 * (this->fvbar.N);
	long beta = inverse(alpha, twoN);
	if (1 != (beta * alpha % twoN)){
		cout << "ERROR: beta=" << beta << " is not the inverse of alpha=" 
			 << alpha << " mod " << twoN << endl;
		exit(11);
	}
	gen_format_fixing_key(this->ffk_a_x, this->ffk_b_x, fvx, fvbar, log_Px, log_qx, beta, crt_fft);
	gen_format_fixing_key(this->ffk_a_y, this->ffk_b_y, fvy, fvbar, log_Py, log_qy, 1, crt_fft);


	cout << "finished gen_format_fixing_keys" << endl;
}

fmpz_polyxx HomomorphicHeatmap::mult_powers_X_by_vec(
									const std::vector<long>& u, 
								  	const fmpz_polyxx& b,
									long N,
									long img_size) {
    fmpz_polyxx res(img_size); // allocate space for img_size coefficients
    for(long i = 0; i < N; i++){
		res.set_coeff(u[i], res.get_coeff(u[i]) + b.get_coeff(i));
	}
    return res;
}

void HomomorphicHeatmap::mult_powers_X_by_circ_mat_by_diag(
												vector<fmpz_polyxx>& res,
												const std::vector<long>& u, 
											  	const fmpz_polyxx& a,
												long img_size,
												const fmpzxx& q,
												char x_or_y
											  	){

	int N;
	vector<vector<fmpzxx> >& res_vec = (x_or_y == 'x' ? res_vec_x : res_vec_y);
	if ('x' == x_or_y){
		N = this->Nx;
	}else{
		N = this->Ny;
	}
	for (int i = 0; i < img_size; i++){
		for (int j = 0; j < N; j++){
			res_vec[i][j] = fmpzxx(0);
		}
	}

	for(int i = 0; i < N; i++){

        const fmpzxx a_i(a.get_coeff(i));

        // add the positive coefficients
		for (int j = i; j < N; j++){
            // x^u[j] * a_i
            res_vec[u[j]][j-i] += a_i;
		}

        // add the negative (rotated) coefficients
		for (int j = 0; j < i; j++){
            // x^u[j] * a_i
//			res_vec[N-i+j][u[j]] -= a_i;
			res_vec[u[j]][N-i+j] -= a_i;
		}
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < img_size; j++){
			res[i].set_coeff(j, res_vec[j][i] % q);
		}
	}

}

MixedCiphertext HomomorphicHeatmap::mult_by_test_vector(const Ciphertext& c,
									const std::vector<long>& func,
									long N,
									long img_size,
									const fmpzxx& q,
									char x_or_y) {
	fmpz_polyxx b = mult_powers_X_by_vec(func, c.b, N, img_size);
	std::vector<fmpz_polyxx> a(N);
	mult_powers_X_by_circ_mat_by_diag(a, func, c.a, img_size, q, x_or_y);
	return MixedCiphertext(a, b);
}


MixedCiphertext HomomorphicHeatmap::apply_func_to_exponent(
											const std::vector<Ciphertext>& cs,
											char x_or_y){
	long k, N, img_size;
	fmpzxx q;

	if ('x' == x_or_y){
		N = this->Nx;
		k = this->k_x;
		q = this->qx;
		img_size = this->img_size_f; //ceil(this->x_max / base_cell);
	}else if ('y' == x_or_y){
		N = this->Ny;
		k = this->k_y;
		q = this->qy;
		img_size = this->img_size_g;//ceil(this->y_max / height_cell);
	}else {
		cout << "ERROR: argument char x_or_y must be 'x' or 'y'." << endl;
		exit(12);
	}


	const std::vector<long>& u_f = ('x' == x_or_y ? this->func_f : this->func_g);
	vector<fmpz_polyxx> a(N);
    fmpz_polyxx b;
	vector<long> u_f_i(N); /* XXX: no need to copy to this vector. We must a parameter
							  		index to the other functions so that they can
									receive this (k*N)-dimensional vector func
									and use only the subvector from index*N to
									(index + 1)*N-1. */

    for (long i = 0; i < k; i++){
		// XXX: no need to copy this. See comment above.
        //u_f_i = func[i*N:(i+1)*N]
		for (long j = 0; j < N; j++)
        	u_f_i[j] = u_f[i*N + j];

        MixedCiphertext u_c_i = mult_by_test_vector(cs[i], u_f_i, N, img_size, q, x_or_y);
        // a += u_c_i.a; // add vectors entrywise
		for (long j = 0; j < N; j++)
			a[j] += u_c_i.a[j];
        b += u_c_i.b;
	}

    for (long i = 0; i < N; i++){
		reduce_poly_mod(a[i], q);
	}
	reduce_poly_mod(b, q);

    return MixedCiphertext(a, b);
}

Ciphertext HomomorphicHeatmap::apply_format_fixing_key(const MixedCiphertext& c, char x_or_y){

	fmpzxx P;
	fmpzxx q;
	std::vector<PolyCRT_FFT> & ffk_a = ('x' == x_or_y ? this->ffk_a_x : this->ffk_a_y);
	std::vector<PolyCRT_FFT> & ffk_b = ('x' == x_or_y ? this->ffk_b_x : this->ffk_b_y);

	long N;

	if ('x' == x_or_y){
		P = this->Px;
		q = this->qx;
		N = this->Nx;
	}else{
		if ('y' == x_or_y){
			P = this->Py;
			q = this->qy;
			N = this->Ny;
		}else{
			cout << "ERROR: Trying to apply inexistent format fixing key." << endl;
			exit(8);
		}
	}

    fmpz_polyxx _b(P * c.b);

	vector<PolyCRT_FFT> fft_c_a(N);
	assert(c.a.size() == N);
	for (long i = 0; i < N; i++){
		fft_c_a[i] = crt_fft->forward(polyFlint_to_Poly64(c.a[i]));
	}
	// inner product: c.a times the "b term" of ffk
	Poly64 _b_k_ = crt_fft->inner_product_full_precomp<int64_t>(fft_c_a, ffk_b);
	fmpz_polyxx _b_k = Poly64_to_PolyFlint(_b_k_);

    fmpz_polyxx b_P(_b - _b_k);

	// inner product: -c.a times the "a term" of ffk
	Poly64 _a_P_ = crt_fft->inner_product_full_precomp<int64_t>(fft_c_a, ffk_a);
	fmpz_polyxx a_P = Poly64_to_PolyFlint(_a_P_);

    // Now, (a_P, b_P) is an encryption of x^f(m) scaled by P and modulo P*q
    // so we modulus-switch it to q by dividing by P
	fmpz_polyxx a_final = rounded_division(a_P, P);
	fmpz_polyxx b_final = rounded_division(b_P, P);

    return Ciphertext(a_final, b_final);
}

	/* cs_x: k ciphertexts as generated by enc_high_level.
	 * cs_y: k ciphertexts as generated by enc_high_level.
	 * Assuming that (cs_x, cs_y) encrypts a point (x, y), this function
	 * returns an RLWE encryption of X^(f(x) * alpha + g(y))
	 */
NonRelinCiphertext HomomorphicHeatmap::compute_cell(
											const vector_ciphertext & cs_x,
											const vector_ciphertext & cs_y
											){

#if VERBOSE
	cout << "applying func to exponent" << endl;
	_t_before = high_resolution_clock::now();
#endif

	MixedCiphertext mxd_c = apply_func_to_exponent(cs_x, 'x');

#if VERBOSE
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "    time to compute enc(X^f(x)): " << _t_duration / 1000.0 << " s" << endl << endl;
#endif


	cout << "     fixing format of x coordinate..." << endl;
#if VERBOSE
	_t_before = high_resolution_clock::now();
#endif

	Ciphertext c_f = apply_format_fixing_key(mxd_c, 'x'); // == enc(X^f(x)) in Rbar under fvbar.sk(X^beta)

#if VERBOSE
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "         time fix format: " << _t_duration / 1000.0 << " s" << endl << endl;
#endif



#if VERBOSE
	_t_before = high_resolution_clock::now();
#endif

	Ciphertext c_f_alpha = fvbar.apply_automorphism(c_f, false); // == enc(X^(alpha * f(x))) under fvbar.sk
	
#if VERBOSE
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "         time to apply automorphism X --> X^alpha: " << _t_duration / 1000.0 << " s" << endl << endl;
#endif



	//  Now, do almost the same for the y coordinate

#if VERBOSE
	_t_before = high_resolution_clock::now();
#endif

	mxd_c = apply_func_to_exponent(cs_y, 'y');

#if VERBOSE
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "    time to compute enc(X^g(y)): " << _t_duration / 1000.0 << " s" << endl << endl;
#endif




	cout << "     fixing format of y coordinate..." << endl;

#if VERBOSE
	_t_before = high_resolution_clock::now();
#endif

	Ciphertext c_g = apply_format_fixing_key(mxd_c, 'y'); // == enc(X^g(y)) in Rbar

#if VERBOSE
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "         time fix format: " << _t_duration / 1000.0 << " s" << endl << endl;
#endif



#if VERBOSE
	_t_before = high_resolution_clock::now();
//	Ciphertext c_f_g = fvbar.mult(c_f_alpha, c_g); // enc(X^(f(x) * alpha + g(y)))
	NonRelinCiphertext c_f_g = fvbar.mult_without_relin(c_f_alpha, c_g); // enc(X^(f(x) * alpha + g(y)))
	_t_after = high_resolution_clock::now();
	_t_duration = duration_cast<milliseconds>(_t_after - _t_before).count();
	cout << "    time to multiply x^f(m) by x^g(m): " << _t_duration / 1000.0 << " s" << endl << endl;
	return c_f_g;
#else
	return fvbar.mult_without_relin(c_f_alpha, c_g); // enc(X^(f(x) * alpha + g(y)))
#endif
}
	
NonRelinCiphertext HomomorphicHeatmap::compute_heatmap(
										const std::vector< vector_ciphertext >& cs_x,
										const std::vector< vector_ciphertext >& cs_y
											){
	if (cs_y.size() != cs_x.size()){
		cout << "ERROR: Trying to compute heatmap with different number of x and y coordinates." << endl;
		exit(11);
	}
	cout << "coordinate (x_0, y_0)..." << endl;
	NonRelinCiphertext heatmap = compute_cell(cs_x[0], cs_y[0]);
	for(long i = 1; i < cs_x.size(); i++){
		cout << "coordinate (x_" << i << ", y_" << i << ")..." << endl;
		NonRelinCiphertext tmp = compute_cell(cs_x[i], cs_y[i]);
		heatmap.a += tmp.a;
		heatmap.b += tmp.b;
		heatmap.c += tmp.c;
	}
	reduce_poly_mod(heatmap.a, fvbar.q);
	reduce_poly_mod(heatmap.b, fvbar.q);
	reduce_poly_mod(heatmap.c, fvbar.q);
	return heatmap;
}


std::ostream& operator<< (std::ostream &out, const HomomorphicHeatmap& hhm) {  
	out << "HomomorphicHeatmap: {" 
		<< "x_max: " << hhm.x_max
		<< ", y_max: " << hhm.y_max
		<< ", base_cell: " << hhm.base_cell
		<< ", height_cell: " << hhm.height_cell
		<< ", img_size_f: " << hhm.img_size_f
		<< ", img_size_g: " << hhm.img_size_g
		<< ", final_img_size: " << hhm.final_img_size
		<< ", alpha: " << hhm.alpha
	    << ", Nx: 2^" << hhm.log_Nx
	    << ", Ny: 2^" << hhm.log_Ny
	    << ", Nbar: 2^" << hhm.log_Nbar
		<< ", k_x: " << hhm.k_x
		<< ", k_y: " << hhm.k_y
	    << ", q_x: 2^" << hhm.log_qx
	    << ", q_y: 2^" << hhm.log_qy
	    << ", q_bar: 2^" << hhm.fvbar.logq
	   << ", t: " << hhm.t
	   << "}";
	return out;
}


#include "HomomorphicHeatmap.h"
#include <cstring>

/* n = 2^9 = 512, log2(q) = 32, sigma = 3000.000000, uniform sk
 * usvp: rop: ~2^128.1,  red: ~2^128.1,  δ_0: 1.004488,  β:  336,  d: 1557,  m:     1044
 * dec: rop: ~2^145.7,  m:     1127,  red: ~^145.7,  δ_0: 1.004170,  β:  375,  d: 1639,  babai~2^131.4,  babai_op: ~2^146.5,  repeat:       72,  ε: 0.062500
 * dual: rop: ~2^140.5,  m:     1142,  red: ~2^140.5,  δ_0: 1.004158,  β:  376,  repeat:~2^97.0,  d: 1654,  c: 
 */


/* n = 2^12 = 4096, log2(q) = 64, sigma = 3.200000, ternary sk
* usvp: rop: ~2^232.9,  red: ~2^232.9,  δ_0: 1.002705,  β:  687,  d: 7843,  m~2^11.9,  repeat:  1,  k:        0,  postprocess:        0
* dec: rop: ~2^258.3,  m: ~2^12.1,  red: ~2^258.3,  δ_0: 1.002501,  β:  764,  d: 8436,  babai~2^249.6,  babai_op: ~2^264.7,  repeat: 7,  ε: 0.500000
* dual: rop: ~2^239.4,  m: ~2^12.0,  red: ~2^239.4,  δ_0: 1.002644,  β:  709,  repeat~2^108.3, d: 8102,  c:  3.919,  k: 64,  postproces  30
*/


using namespace std;

using namespace std::chrono; // to measure execution times

vector<long> random_vector(int dimension, int mod){
	vector<long> f(dimension);
	for(int i = 0; i < f.size(); i++){
		f[i] = random() % mod; //RandomBnd(img_size);
	}
	return f;
}



vector<long> random_func(int domain_size, int img_size){
	return random_vector(domain_size, img_size);
}


void test_split_domain_heatmap(long x_max, long base_cell) {

	// To measure run time
	high_resolution_clock::time_point t_before;
    high_resolution_clock::time_point t_after;
	double t_duration;
	
	// Variables of the heapmap
	long y_max = x_max;
	long height_cell = base_cell;

	long log_Nx = 9;
	long log_Ny = 9;

	long log_qx = 32;
	long log_qy = 32;

	long t = (1<<2);

	int npoints = 10;

	// defining input points. 
	vector<long> x_coords = random_vector(npoints, x_max);
	vector<long> y_coords = random_vector(npoints, y_max);


	cout << "Generating keys... " << endl;
	t_before = high_resolution_clock::now();
	HomomorphicHeatmap hhm = HomomorphicHeatmap(x_max, y_max, base_cell, height_cell, log_Nx, log_Ny, log_qx, log_qy, t);
	t_after = high_resolution_clock::now();
	t_duration = duration_cast<milliseconds>(t_after - t_before).count();
	cout << "    time to generate keys: " << t_duration / 1000.0 << " s" << endl << endl;

	cout << hhm << endl;

	cout << "Encrypting x_coords ..." << endl;
	t_before = high_resolution_clock::now();
	vector<vector_ciphertext > cs_x = hhm.enc_high_level(x_coords, 'x');
	t_after = high_resolution_clock::now();
	t_duration = duration_cast<milliseconds>(t_after - t_before).count();
	cout << "    time to encrypt: " << t_duration << " ms" << endl << endl;
	
	cout << "Encrypting y coords ..." << endl;
	t_before = high_resolution_clock::now();
	vector<vector_ciphertext > cs_y = hhm.enc_high_level(y_coords, 'y');
	t_after = high_resolution_clock::now();
	t_duration = duration_cast<milliseconds>(t_after - t_before).count();
	cout << "    time to encrypt: " << t_duration << " ms" << endl << endl;


	cout << "Computing heapmap ..." << endl;
	t_before = high_resolution_clock::now();

	NonRelinCiphertext c = hhm.compute_heatmap(cs_x, cs_y);

	t_after = high_resolution_clock::now();
	t_duration = duration_cast<milliseconds>(t_after - t_before).count();
	cout << "Time to compute heatmap: " << t_duration/1000.0 << " s"  << endl << endl;
	cout << "Time per point: " << t_duration/(x_coords.size() * 1000.0) << " s"  << endl << endl;

	fmpz_polyxx dec_msg = hhm.fvbar.dec(c);


	string temp = "x";
	char x[2];
	strcpy(x, temp.c_str());
//	cout << dec_msg << endl;
	cout << "Computed heatmap: " << dec_msg.pretty(x) << endl;

	fmpz_polyxx expected;
	for (long i = 0; i < x_coords.size(); i++){
		long j = floor(x_coords[i] / (double)base_cell) * hhm.alpha + floor(y_coords[i] / (double) height_cell);
		expected.set_coeff(j, expected.get_coeff(j) + 1);
	}

	cout << "Final noise: " << hhm.fvbar.get_noise(c, dec_msg) << endl;

	if (expected == dec_msg){
		cout << "test_homomorphic_heatmap: OK" << endl;
        cout << "----------------------------------" << endl;
        cout << endl;
	}else{
		cout << "test_homomorphic_heatmap: NOT OK" << endl;
		cout << "computed: " << dec_msg.pretty(x) << endl;
		cout << "expected: " << expected.pretty(x) << endl;
        exit(1);
	}
}

int main() {

    test_split_domain_heatmap(1 << 10, 1 << 6);
    test_split_domain_heatmap(1 << 11, 1 << 7);
    test_split_domain_heatmap(1 << 12, 1 << 8);
    test_split_domain_heatmap(1 << 13, 1 << 9);
    test_split_domain_heatmap(1 << 14, 1 << 9);
    test_split_domain_heatmap(1 << 15, 1 << 9);

	return 0;
}

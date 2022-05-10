#include "testvectors.h"
#include <bitset>
#include <string>
#include <sstream>
#include <chrono>


using namespace std;
using namespace seal;

vector<long> random_vector(int dimension, int mod){
	vector<long> v(dimension);
    for(int i = 0; i < v.size(); i++){
        v[i] = random() % mod;
    }
	return v;
}


int func(int x, char xy, int base_cell, int height_cell, int max_y) {
    if (xy == 'x') {
        return floor(x/base_cell)*ceil(max_y/height_cell);
    } else if (xy == 'y') {
        return floor(x/height_cell);
    } else {
        throw std::invalid_argument("xy should be either \'x\' or \'y\'.");
    }
}

void get_x_power_f_m(seal::Ciphertext ct_m, int nb_bits, Plaintext T_plain[], seal::Ciphertext &dest, seal::EncryptionParameters parms, seal::GaloisKeys galois_keys, seal::RelinKeys relin_keys, int N_inverse, SEALContext& context, Evaluator& evaluator) {

    const int t = *data(parms.plain_modulus());
    const int N = parms.poly_modulus_degree();
    
    std::ostringstream N_stream;
    N_stream << std::hex << N_inverse;
    Plaintext N_inv(N_stream.str());
    std::ostringstream str_stream;
    str_stream << std::hex << t-1;
    std::string n_one_hex =  str_stream.str();    
    
    //cout << "Multiplying testvectors with x^z..." << endl;
    Ciphertext x_fx[nb_bits];
    auto start = std::chrono::high_resolution_clock::now();
    Plaintext zero("0");
    int fixed_bits[nb_bits]; // indicates when a bit is always zero.
    for (int i=0; i<nb_bits; i++) {
        Plaintext& T_plain_i = T_plain[i];
        if (T_plain_i == zero) {
            fixed_bits[i] = 1;
        } else {
            evaluator.multiply_plain(ct_m, T_plain_i, x_fx[i]);
            fixed_bits[i] = 0;
        }
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time to mulitply by X^z: " << elapsed.count() << endl;
   
 
    //cout << "Automorphisms to extract constant term..." << endl;
    start = std::chrono::high_resolution_clock::now();
    int logN = ceil(log2(N));
    Ciphertext temp;
    for (int i = 0; i<nb_bits; i++) {
        if (fixed_bits[i] == 0) {
            for (int j = 0; j<logN; j++) {
                evaluator.apply_galois(
                    x_fx[i],
                    (N/pow(2,j)+1),
                    galois_keys,
                    temp);
                evaluator.add_inplace(x_fx[i],temp);
            }
            evaluator.multiply_plain_inplace(x_fx[i],N_inv);
            std::string poly = "1x^" + to_string(int(pow(2,nb_bits -1 - i))) + " + " + n_one_hex;
            Plaintext plain_poly(poly);
            Plaintext plain_one("1");
            evaluator.multiply_plain_inplace(x_fx[i],plain_poly);
            evaluator.add_plain_inplace(x_fx[i],plain_one);
        }                    
    } 
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "Elapsed time to apply automorphism to extract constant term: " << elapsed.count() << endl;


    // Now, perform multiplication to group the bits in the exponent of X.
    // We use a tree structure to reduce noise.
    start = std::chrono::high_resolution_clock::now();
    int tree_depth = ceil(log2(nb_bits));
    for (int i=1; i<=tree_depth; i++) {
        for (int j=0; j<nb_bits; j += pow(2,i)) {
            int j2 = j+pow(2,i-1);
            if (nb_bits > j2) {
                if (fixed_bits[j]==0 && fixed_bits[j2] ==0 ) {
                    evaluator.multiply_inplace(x_fx[j], x_fx[j2]);
                    evaluator.relinearize_inplace(x_fx[j], relin_keys);
                } else if (fixed_bits[j]!=0 && fixed_bits[j2]==0 ) {
                    x_fx[j] = x_fx[j2];
                    fixed_bits[j] = 0;
                } 
            }
        }
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "Elapsed time for tree mulitplication: " << elapsed.count() << endl;
    
    dest = x_fx[0];
}

int inverse_modulus(int a, int b) { 
    int q = 0;
    int r = 1;
    int s1 = 1; 
    int s2 = 0;
    int s3 = 1; 
    int t1 = 0; 
    int t2 = 1;   
    int t3 = 0;  
    while(r > 0){
        q = floor(a/b);
	r = a - q * b;
	s3 = s1 - q * s2;
	t3 = t1 - q * t2;
		
	if(r > 0){
	    a = b;
	    b = r;
	    s1 = s2;
	    s2 = s3;
	    t1 = t2;
	    t2 = t3;
	}
    }
    return s2;
}

std::string dec2binstr(int a, int nb_bits)
{
    std::string binary = "";
    int mask = 1;
    for(int i = 0; i < nb_bits; i++)
    {
        if((mask&a) >= 1)
            binary = "1"+binary;
        else
            binary = "0"+binary;
        mask<<=1;
    }
//    cout<<binary<<endl;
    return binary;
}

void get_testvectors(seal::EncryptionParameters parms, char xy, int nb_bits, int max, Plaintext T_plain[], int base_cell, int height_cell, int max_y) 
{
    std::string T[nb_bits];
    const int t = *data(parms.plain_modulus());
    const int N = parms.poly_modulus_degree();
    int N_inverse = inverse_modulus(N, t);
    if (N_inverse < 0 ) {
        N_inverse += t;
    } 

    std::ostringstream N_stream;
    N_stream << std::hex << N_inverse;
    Plaintext N_inv(N_stream.str());
    std::ostringstream str_stream;
    str_stream << std::hex << t-1;
    std::string n_one_hex =  str_stream.str();
    
    int LUTi = func(1, xy, base_cell, height_cell, max_y);
    std::string bitdec1 = dec2binstr(LUTi,nb_bits);
    for (int j=0; j<nb_bits; j=j+1) {
        if (bitdec1[j] != '0') {
            T[j] = T[j] + n_one_hex + "x^" + to_string(N-1);
        }
        else {
            T[j] = T[j] + "0x^" + to_string(N-1);
        }
    }
    for (int i=2; i<max; i=i+1) {
        LUTi = func(i, xy, base_cell, height_cell, max_y);
        std::string bitdeci = dec2binstr(LUTi,nb_bits);
        for (int j=0; j<nb_bits; j=j+1) {
            if (bitdeci[j] != '0') {
                T[j] = T[j] + " + " + n_one_hex + "x^" + to_string(N-i);
            }
        }
    }
    LUTi = func(0, xy, base_cell, height_cell, max_y);
    std::string bitdec0 = dec2binstr(LUTi,nb_bits);
    for (int j=0; j<nb_bits; j=j+1) {
        if (bitdec0[j] != '0') {
            T[j] = T[j] + " + " + n_one_hex;
        }
    }

	// convert to seal::Plaintext
	for(int i = 0; i < nb_bits; i++)
		T_plain[i] = Plaintext(T[i]);

}


void test_split_domain_heatmap(int max_x, int base_cell)
{

    cout << "max_x = 2^" << log(max_x)/log(2) << endl;
    cout << "base_cell = 2^" << log(base_cell)/log(2) << endl;

    int max_y = max_x;
    int height_cell = base_cell;

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 1<<12;
    if (poly_modulus_degree < max_x)
        poly_modulus_degree = max_x;

    parms.set_poly_modulus_degree(poly_modulus_degree);

    // Now, depending on N, choose the modulus Q with the number of bits decribed
    // in Table 1 of the paper
    if ((1<<12) == poly_modulus_degree)
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    else if ((1<<13) == poly_modulus_degree)
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 50, 40, 40}));
    else if ((1<<14) == poly_modulus_degree)
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 50, 50, 45}));
    else 
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 55, 52, 50}));

	Modulus t(5);

    parms.set_plain_modulus(t);
    SEALContext context(parms);
    print_parameters(context);
    
	   
	int npoints = 10;
    /*  INPUTS ARE DEFINED HERE  */
	vector<long> x = random_vector(npoints, max_x);
	vector<long> y = random_vector(npoints, max_y);

    if (max_x > poly_modulus_degree || max_y > poly_modulus_degree) {
        throw std::invalid_argument("The parameters for this heatmap are incorrect. The degree N must be larger than x_max and y_max.");
    }
 
    int ncells = floor((max_x-1)/base_cell)*ceil((max_y)/height_cell) + floor((max_y-1) / height_cell);
	cout << "Number of cells: " << ncells << " = 2^" << log(ncells) / log(2) << endl;
    if (ncells > poly_modulus_degree) {
        throw std::invalid_argument("The parameters for this heatmap are incorrect. The degree N must be larger than the number of cells.");
    }
    
    cout << "Parameter validation (success): " << context.parameter_error_message() << endl;
    
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    vector<uint32_t> gal_elts;
    int logN = ceil(log2(poly_modulus_degree));
    for (int i=0; i<logN; i++) {
        gal_elts.push_back((uint32_t) (poly_modulus_degree/pow(2,i)+1));
    }
    GaloisKeys galois_keys;
    galois_keys = keygen.create_galois_keys(gal_elts,false);
    
    

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);


    const int t_int = *data(parms.plain_modulus());
    const int N = parms.poly_modulus_degree();
    int N_inverse = inverse_modulus(N, t_int);
    if (N_inverse < 0 ) {
        N_inverse += t_int;
    }

    Ciphertext heatmap;
    Ciphertext x_coord;
    Ciphertext y_coord;
    Plaintext pt_x("1x^"+to_string(x[0]));
    Plaintext pt_y("1x^"+to_string(y[0]));
    Ciphertext ct_x;
    Ciphertext ct_y;
    char char_x = 'x';
    char char_y = 'y'; 
    const int nb_bits_x = ceil(log2(func(max_x,char_x, base_cell, height_cell, max_y)));
    const int nb_bits_y = ceil(log2(func(max_y,char_y, base_cell, height_cell, max_y)));

    std::string T_x[nb_bits_x], T_y[nb_bits_y];
	Plaintext T_plain_x[nb_bits_x];
	Plaintext T_plain_y[nb_bits_y];
    get_testvectors(parms,char_x,nb_bits_x,max_x,T_plain_x, base_cell, height_cell, max_y);
    get_testvectors(parms,char_y,nb_bits_y,max_y,T_plain_y, base_cell, height_cell, max_y);

    encryptor.encrypt(pt_x, ct_x); 
    encryptor.encrypt(pt_y, ct_y);
    get_x_power_f_m(ct_x, nb_bits_x, T_plain_x, x_coord, parms, galois_keys, relin_keys,N_inverse, context, evaluator);
    get_x_power_f_m(ct_y, nb_bits_y, T_plain_y, y_coord, parms, galois_keys, relin_keys,N_inverse, context, evaluator);
    evaluator.multiply(x_coord, y_coord, heatmap);
    
	double avg_time_per_point = 0;

    for (int i = 1; i< x.size(); i++) {
        cout << "Processing point " << i  << endl;
        
        Plaintext pt_x("1x^"+to_string(x[i]));
        Plaintext pt_y("1x^"+to_string(y[i]));
        Ciphertext ct_x;
        Ciphertext ct_y;
        encryptor.encrypt(pt_x, ct_x); 
        encryptor.encrypt(pt_y, ct_y); 
       
        auto start_tot = std::chrono::high_resolution_clock::now();
        Ciphertext res_ct;


        get_x_power_f_m(ct_x, nb_bits_x, T_plain_x, ct_x, parms, galois_keys, relin_keys,N_inverse, context, evaluator);
        get_x_power_f_m(ct_y, nb_bits_y, T_plain_y, ct_y, parms, galois_keys, relin_keys,N_inverse, context, evaluator);

        auto start_mult = std::chrono::high_resolution_clock::now();

        evaluator.multiply_inplace(ct_x, ct_y);

        auto finish_mult = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_mult = finish_mult - start_mult;
		double _elapsed_mult = elapsed_mult.count();
        cout << "Elapsed time for hom. mult: " << _elapsed_mult << endl;

        evaluator.add_inplace(heatmap, ct_x);
        
        auto finish_tot = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish_tot - start_tot;
		double elapsed_point = elapsed.count();
        cout << "Total elapsed time for this point: " << elapsed_point << endl;
		cout << endl;

		avg_time_per_point += elapsed_point;
    }
	avg_time_per_point /= (npoints - 1);
    

    Plaintext res;
    decryptor.decrypt(heatmap,res);
    double noise_budget = decryptor.invariant_noise_budget(heatmap);
    
    // a minimal acceptable noise bound
    if (noise_budget <= 2){
        cout << "ERROR: too much noise in the final result." << endl;
        cout << "Noise budget: " << noise_budget << " bits" << endl;
        exit(1);
    }


    cout << "result: " << res.to_string() << endl;
    cout << "noise budget in encrypted_result: " << noise_budget << " bits" << endl;
	cout << "AVG time per point: " << avg_time_per_point << endl;
    cout << "------------------------------" << endl;
    cout << endl;
}



int main(){
    test_split_domain_heatmap(1 << 10, 1 << 6);
    test_split_domain_heatmap(1 << 11, 1 << 7);
    test_split_domain_heatmap(1 << 12, 1 << 8);
    test_split_domain_heatmap(1 << 13, 1 << 9);
    test_split_domain_heatmap(1 << 14, 1 << 9);
    test_split_domain_heatmap(1 << 15, 1 << 9);
    return 0;
}

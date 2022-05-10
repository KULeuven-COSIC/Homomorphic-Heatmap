#ifndef __FFT_engine_
#define __FFT_engine_

#include <vector>
#include <fftw3.h>
#include <complex>
#include <map>

using namespace std;

typedef std::vector<int32_t> Poly32;
typedef std::vector<int64_t> Poly64;
typedef std::vector<double> PolyDouble;
typedef std::vector<std::complex<double>> PolyFFT;

//TODO: description

class FFT_engine
{
    int fft_dim;
    int fft_dim2;

    fftw_plan plan_to_fft;
    fftw_plan plan_from_fft;

public:
    double* in_array;
    fftw_complex* out_array;
    
	void from_fft_core(const PolyFFT& in) const;

    vector<PolyFFT> pos_powers;
    vector<PolyFFT> neg_powers;

    FFT_engine() = delete;
    FFT_engine(const int dim);

    void to_fft(PolyFFT& out, const Poly32& in) const;
    void to_fft(PolyFFT& out, const Poly64& in) const;
    void to_fft(PolyFFT& out, const PolyDouble& in) const;
    void from_fft(Poly32& out, const PolyFFT& in) const;
    void from_fft(Poly64& out, const PolyFFT& in) const;
    void from_fft(PolyDouble& out, const PolyFFT& in) const;

    ~FFT_engine();
};

PolyFFT operator *(const PolyFFT& a, const PolyFFT& b);
void operator *=(PolyFFT& a, const PolyFFT& b);
PolyFFT operator *(const PolyFFT& a, const int b);
PolyFFT operator +(const PolyFFT& a, const PolyFFT& b);
void operator +=(PolyFFT& a, const PolyFFT& b);
void operator +=(PolyFFT& a, const complex<double> b);
PolyFFT operator -(const PolyFFT& a, const PolyFFT& b);
void operator -=(PolyFFT& a, const PolyFFT& b);

void operator +=(Poly32& a, const Poly32& b);
void operator +=(Poly64& a, const Poly64& b);
#endif

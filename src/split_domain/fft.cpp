#include "fft.h"
#include <cassert>
#include <iostream>

FFT_engine::FFT_engine(const int dim): fft_dim(dim)
{
    assert(dim%2 == 0);

    fft_dim2 = (dim >> 1) + 1;

    in_array = (double*) fftw_malloc(sizeof(double) * 2*dim);
    out_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (dim + 2));
    plan_to_fft = fftw_plan_dft_r2c_1d(2*dim, in_array, out_array,  FFTW_PATIENT);
    plan_from_fft = fftw_plan_dft_c2r_1d(2*dim, out_array, in_array,  FFTW_PATIENT);

    pos_powers = vector<PolyFFT>(dim,PolyFFT(fft_dim2));
    neg_powers = vector<PolyFFT>(dim,PolyFFT(fft_dim2));
    for(int i = 0; i < dim; i++)
    {
        Poly32 x_power(dim,0);
        //x_power[0] = -1;
        x_power[i] += 1;
        PolyFFT x_power_fft;
        to_fft(x_power_fft, x_power);
        pos_powers[i] = x_power_fft;

        x_power[i] -= 2;
        to_fft(x_power_fft, x_power);
        neg_powers[i] = x_power_fft;
    }
    //x_powers.insert({{-1,neg_powers}, {1,pos_powers}});
}

void FFT_engine::to_fft(PolyFFT& out, const Poly32& in) const
{
    for (int i = 0; i < fft_dim; i++)
    {
        in_array[i] = double(in[i]);
        in_array[i+fft_dim] = 0.0;
    }
    fftw_execute(plan_to_fft);
    out.resize(fft_dim2);
    int tmp = 1;
    for (int i = 0; i < fft_dim2; i++)
    {
        fftw_complex& out_arr = out_array[tmp];
        complex<double>& outi = out[i];
        outi.real(out_arr[0]);
        outi.imag(out_arr[1]);
        tmp += 2;
    }
}

void FFT_engine::to_fft(PolyFFT& out, const Poly64& in) const
{
    for (int i = 0; i < fft_dim; i++)
    {
        in_array[i] = double(in[i]);
        in_array[i+fft_dim] = 0.0;
    }
    fftw_execute(plan_to_fft);
    out.resize(fft_dim2);
    int tmp = 1;
    for (int i = 0; i < fft_dim2; i++)
    {
        fftw_complex& out_arr = out_array[tmp];
        complex<double>& outi = out[i];
        outi.real(out_arr[0]);
        outi.imag(out_arr[1]);
        tmp += 2;
    }
}

void FFT_engine::to_fft(PolyFFT& out, const PolyDouble& in) const
{
    for (int i = 0; i < fft_dim; i++)
    {
        in_array[i] = double(in[i]);
        in_array[i+fft_dim] = 0.0;
    }
    fftw_execute(plan_to_fft);
    out.resize(fft_dim2);
    int tmp = 1;
    for (int i = 0; i < fft_dim2; i++)
    {
        fftw_complex& out_arr = out_array[tmp];
        complex<double>& outi = out[i];
        outi.real(out_arr[0]);
        outi.imag(out_arr[1]);
        tmp += 2;
    }
}


void FFT_engine::from_fft_core(const PolyFFT& in) const
{
    int tmp = 0;
    for (int i = 0; i < fft_dim2; i++) 
    {
        //std::cout << "i: " << i << ", number: " << in[i] << std::endl;
        out_array[tmp+1][0] = real(in[i])/fft_dim;
        out_array[tmp+1][1] = imag(in[i])/fft_dim;
        out_array[tmp][0] = 0.0;
        out_array[tmp][1] = 0.0;
        tmp += 2;
    }
    fftw_execute(plan_from_fft);
}

void FFT_engine::from_fft(PolyDouble& out, const PolyFFT& in) const
{
	from_fft_core(in);

    out.resize(fft_dim);
    for (int i = 0; i < fft_dim; i++)
    {	
        out[i] = in_array[i];
    }
}


void FFT_engine::from_fft(Poly64& out, const PolyFFT& in) const
{
	from_fft_core(in);


    out.resize(fft_dim);
    for (int i = 0; i < fft_dim; i++)
    {	
        out[i] = (int64_t)(round(in_array[i]));
        //std::cout << "i: " << i << ", number: " << out[i] << std::endl;
    }
}


void FFT_engine::from_fft(Poly32& out, const PolyFFT& in) const
{
	from_fft_core(in);

    out.resize(fft_dim);
    for (int i = 0; i < fft_dim; i++)
    {	
        out[i] = (int32_t)(round(in_array[i]));
        //std::cout << "i: " << i << ", number: " << out[i] << std::endl;
    }
}

FFT_engine::~FFT_engine()
{
    fftw_destroy_plan(plan_to_fft);
    fftw_destroy_plan(plan_from_fft);
    fftw_free(in_array);
    fftw_free(out_array);
}

PolyFFT operator +(const PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    PolyFFT res(a.size());
    for (size_t i = 0; i < a.size(); i++)
        res[i] = a[i]+b[i];

    return res;
}

void operator +=(PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i]+=b[i];
}

void operator +=(PolyFFT& a, const complex<double> b)
{
    for (size_t i = 0; i < a.size(); i++)
        a[i]+=b;
}

PolyFFT operator -(const PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    PolyFFT res(a.size());
    for (size_t i = 0; i < a.size(); i++)
        res[i] = a[i]-b[i];

    return res;
}

void operator -=(PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i]-=b[i];
}

PolyFFT operator *(const PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    PolyFFT res(a.size());
    for (size_t i = 0; i < a.size(); i++)
        res[i] = a[i]*b[i];

    return res;
}

void operator *=(PolyFFT& a, const PolyFFT& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i]*=b[i];
}

// TODO: make a test
PolyFFT operator *(const PolyFFT& a, const int b)
{
    PolyFFT res(a.size());
    double bd = double(b);
    for (size_t i = 0; i < a.size(); i++)
        res[i] = a[i] * bd;

    return res;
}



void operator +=(Poly32& a, const Poly32& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i]+=b[i];
}

void operator +=(Poly64& a, const Poly64& b)
{
    // check that input vectors have the same size
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i]+=b[i];
}

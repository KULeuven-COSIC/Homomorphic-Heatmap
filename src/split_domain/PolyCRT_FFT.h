#ifndef __POLY_CRT_FFT__
#define __POLY_CRT_FFT__

#include "fft.h"

typedef std::vector<PolyFFT > PolyCRT_FFT;

/* Operations on PolyCRT_FFT */
void operator +=(PolyCRT_FFT& a, const PolyCRT_FFT& b);
	
void operator *=(PolyCRT_FFT& a, const PolyCRT_FFT& b);


#endif

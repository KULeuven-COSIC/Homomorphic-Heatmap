
#include "PolyCRT_FFT.h"
#include <cassert>

/* Operations on PolyCRT_FFT */
void operator +=(PolyCRT_FFT& a, const PolyCRT_FFT& b) {
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i] += b[i];
}
void operator *=(PolyCRT_FFT& a, const PolyCRT_FFT& b) {
    assert(a.size() == b.size());

    for (size_t i = 0; i < a.size(); i++)
        a[i] *= b[i];
}



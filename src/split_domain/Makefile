NTL_FLAGS=-I/usr/local/include/  -L/usr/local/lib -lntl
GMP_FLAGS=-lgmp
FFTW_FLAGS=-lfftw3
FLINT_FLAGS=-I/usr/local/include/ -L/usr/local/lib -lflint

CCXFLAGS=-O3 -funroll-loops -march=native -std=c++17 -pthread -I. -I./include -I/usr/local/include/
LIBS=-lflint -lntl -lgmp -lfftw3 -lm

all: test_FV test_HomomorphicHeatmap

utils.o: utils.h utils.cpp
	g++ utils.h utils.cpp -c -g $(CCXFLAGS) $(NTL_FLAGS) $(GMP_FLAGS) $(FLINT_FLAGS) $(LIBS)

FV.o: FV.h FV.cpp utils.o
	g++ FV.cpp FV.h -c -g -O3 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) $(FLINT_FLAGS) -lm

HomomorphicHeatmap.o: HomomorphicHeatmap.h HomomorphicHeatmap.cpp FV.o utils.o crt_fft.o
	g++ HomomorphicHeatmap.cpp HomomorphicHeatmap.h -c -g $(CCXFLAGS) $(NTL_FLAGS) $(GMP_FLAGS) $(FLINT_FLAGS) -lm

test_HomomorphicHeatmap: test_HomomorphicHeatmap.cpp HomomorphicHeatmap.o FV.o utils.o crt_fft.o fft.o
	g++ test_HomomorphicHeatmap.cpp *.o $(CCXFLAGS) $(NTL_FLAGS) $(GMP_FLAGS) $(FLINT_FLAGS) $(FFTW_FLAGS) -o test_HomomorphicHeatmap

fft.o: fft.h fft.cpp
	g++ $(CCXFLAGS) -c fft.cpp $(LIBS)

polyCRT_FFT.o: PolyCRT_FFT.h PolyCRT_FFT.cpp
	g++ $(CCXFLAGS) -c PolyCRT_FFT.cpp $(LIBS) -o polyCRT_FFT.o

crt_fft.o: CRT_FFT.h CRT_FFT.cpp fft.o utils.o polyCRT_FFT.o
	g++ $(CCXFLAGS) $(FLINT_FLAGS) -c CRT_FFT.cpp $(LIBS) -o crt_fft.o

clean:
	rm *.o *.gch test_HomomorphicHeatmap

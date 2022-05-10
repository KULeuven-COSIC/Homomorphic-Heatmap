# Homomorphically counting elements with the same property

Two methods to fast computation of heatmaps over encrypted data using fully 
homomorphic encryption (FHE), as described in the paper
*Homomorphically counting elements with the same property*,
by [Ilia Iliashenko](https://homes.esat.kuleuven.be/~ilia/), 
Malika Izabachène, 
[Axel Mertens](https://www.esat.kuleuven.be/cosic/people/axel-mertens/),
and [Hilder Vitor Lima Pereira](hilder-vitor.github.io/).

## Dependencies

This codes uses [GMP](https://gmplib.org/), [FFTW](https://www.fftw.org/), [FLINT](https://www.flintlib.org/), 
[NTL](https://libntl.org) and [SEAL](https://github.com/microsoft/SEAL/).

To overcome a limitation of SEAL when choosing the parameters, we had to slightly change 
(remove two checks) the files seal/keygenerator.cpp and seal/keygenerator.h 
The modified files can be found in the directory src/install_scripts/our_seal_files/.
The original files from SEAL must be replaced by ours before SEAL is installed.

## Installation

To simplify the installation, we prepared Bash scripts to download and install all
the dependencies in standard directories.
Hence, one can simply run the following command on a Linux terminal

`cd src/ && ./install_third_party_libs.sh`

One can also change our scripts to install the dependencies in local directories,
removing thus the need for root access. 
However, to compile the code, the path to these libs will have to be indicated.

## Compiling and running the code

To compute homomorphic heatmaps using the *full-domain strategy*, just enter in the directory src/ and run

`./run_full_domain.sh`

To compute homomorphic heatmaps using the *split-domain strategy*, just enter in the directory src/ and run

`./run_split_domain.sh`

These scripts will compile and run the code.

If you install the dependencies in non standard directories, then you have to indicate the path in our Makefiles.

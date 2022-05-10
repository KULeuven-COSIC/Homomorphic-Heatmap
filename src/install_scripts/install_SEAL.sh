#!/bin/bash
mkdir SEAL
cd SEAL
wget https://github.com/microsoft/SEAL/archive/refs/heads/main.zip
unzip main.zip
mv SEAL-main/* .

### Replace keygenerator.h and keygenerator.cpp by our slightly modified versions of the files
cp ../our_seal_files/keygenerator.* native/src/seal/

cmake -S . -B build

sudo cmake --build build

cd ..

#!/bin/bash

cd full_domain
cmake -S . -B build
cmake --build build
cd ..
./full_domain/build/bin/testvectors

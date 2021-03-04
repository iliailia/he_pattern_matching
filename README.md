# Homomorphic pattern matching algorithm

This repository contains the source code of the homomorphic pattern matching algorithm from the paper [Homomorphic string search with constant multiplicative depth](https://eprint.iacr.org/2020/931) by Charlotte Bonte and Ilia Iliashenko.

To build the code, install [HElib](https://github.com/homenc/HElib) and run 

    cmake .

and finally

    make

To run the code, use the following syntax
  
    ./pattern_matching p m q pattern_length experiment_runs wildcard_bool

The following line performs 100 tests of the string search algorithm with the parameter set `7.12*` and pattern length 7
  
    ./pattern_matching 7 21177 320 7 100 1

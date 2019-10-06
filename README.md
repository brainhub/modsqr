# Fast modulo squaring for modulo exponentiation
This project provides fast native code to perform sequential modulo squaring on x86_64. This code can also be used as a basis for modulo exponentiation in cryptographic libraries, given that squaring is a dominant portion of modulo exponentiation, which in turn is a building block of many cryptographic algorithms.

Currently the code performs operations in a 2048-bit integer field. This corresponds to a 2048-bit DH exponentiation, 2048-bit [Verifiable Delay Function](https://eprint.iacr.org/2018/627) calculation, or 4096-bit RSA via [RSA CRT method](https://en.wikipedia.org/wiki/RSA_(cryptosystem)#Using_the_Chinese_remainder_algorithm), to name a few.

The project currently offers 4 tools:
1. A squaring library, consisting of a 3-function API in the form Init/Calculate/Free. This library is suitable for integration into other projects.
2. `sqr_test` tool that provides:
   - Benchmarking of modulo squaring performance 
   - Report on x86_64 CPU features that are important to achieve top performance
   - Estimate of the VDF delay parameter `t`, such that it will take 1 day on to perform `2^t` squares on the test machine
3. `sqr` command-line tool to perform sequential squaring. This is the easiest way to interface with this project.
4. A short example on how to integrate `sqr` with JavaScript

## Table of Contents
 - [Background](#background)
 - [Building](#building)
 - [Installation](#installation)
 - [Usage](#usage)
 - [API](#api)
 - [Style](#style)
 - [License](#license)
 
## Background

This project is helpful in _proving_ a [Simple Verifiable Delay Functions](https://eprint.iacr.org/2018/627) value and brute-force unlocking of a secret value. 

Because squaring is the main operation in modulo exponentiation, the code in this project is a good benchmarking tool to time modulo exponentiation operations, such as RSA decryption or encryption, DH in the modulo prime number field, and many more complex protocols that perform exponentiations modulo a large integer.

## Building

After cloning the project and changing the directory into the location of this [README](README.md), run:
```
make
```

## Installation

### Dependencies

This project relies on the low-level OpenSSL library `libcrypto.so`. On Linux Fedora this library can be installed as follows:
```
dnf install openssl-libs
```
The code doesn't depend on the larger OpenSSL package that includes the `openssl` application and other libraries. 

The AVX-512 code path has minimal OpenSSL dependencies and is used for conversion of input and output only. There are no OpenSSL API functions invoked in performance-critical code on the AVX-512 code path.

If you plan to contribute code to this respository, you should install `uncrustify` from [uncrustify](https://github.com/uncrustify/uncrustify) somewhere in your `PATH`.

### Componenets

1. The library `sqrlib.a` with the API interface defined in `sqrlib.h`.
2. `sqr` command line tool to perform production-time calculation
3. `sqr_test` benchmarking and testing tool that can also be invoked as 
   ```
   make test
   ```
4. JavaScript [example](example/vdf-prover.js)

## Usage

### Command-line tool `sqr_test`

Here is a sample output of this tool on an Ivy Bridge CPU architecture: 
```
$ ./sqr_test 
CPU information:
  Brand              : GenuineIntel
  Brand              : Intel(R) Core(TM) i5-3550 CPU @ 3.30GHz
  Features           : ebx=02100800 ecx=7fbae3ff edx=bfebfbff
  Ext. features      : ebx=00000281 ecx=00000000 edx=9c000400
  SSE2               : yes
  AVX                : yes
  MULX, ADX          : no (LOW PERFORMANCE)
  AVX512, AVX512IFMA : no (IMPORTANT!)

OpenSSL version: OpenSSL 1.1.0i-fips  14 Aug 2018

721778 op/sec in 1.45 sec for x^2^2^20
cycles for one square mod N: 4561

Estimated time to complete the next test is 46 seconds
721704 op/sec in 46.49 sec for x^2^2^25
cycles for one square mod N: 4562

x^2^2^36 will take approximately 1 day and 146.972 min (1.10206 day) to compute on this system (t=36)
```

Here is a sample run on a more recent CPU. Notice that the warning with `LOW PERFORMANCE` is gone, but the one with `IMPORTANT!` remains. A higher value of `t` is suggested by the `sqr_test`. 

```
$ make test
./sqr_test
CPU information:
  Brand           : GenuineIntel
  Brand           : Intel(R) Core(TM) i7-6800K CPU @ 3.40GHz
  Features        : ebx=02100800 ecx=7ffefbbf edx=bfebfbff
  Ext. features   : ebx=021cbfbb ecx=00000000 edx=9c000000
  SSE2            : yes
  AVX             : yes
  MULX, ADX       : yes
  AVX512, AVX512IFMA : no (IMPORTANT!)

OpenSSL version: OpenSSL 1.1.1c FIPS  28 May 2019

1.40689e+06 op/sec in 0.75 sec for x^2^2^20
cycles for one square mod N: 2415

Estimated time to complete the next test is 23 seconds
1.50313e+06 op/sec in 22.32 sec for x^2^2^25
cycles for one square mod N: 2260

x^2^2^37 will take approximately 1 day and 83.9173 min (1.05828 day) to compute on this system (t=37)
```

Finally, here is a test on a AVX-512-capable CPU:

```
$ ./sqr_test 
CPU information:
  Brand              : GenuineIntel
  Brand              : Intel(R) Core(TM) i3-8121U CPU @ 2.20GHz
  Features           : ebx=02100800 ecx=7ffafbbf edx=bfebfbff
  Ext. features      : ebx=f2bf67eb ecx=0000001e edx=9c000000
  SSE2               : yes
  AVX                : yes
  MULX, ADX          : yes
  AVX512, AVX512IFMA : yes (BEST PERFORMANCE)

OpenSSL version: OpenSSL 1.1.1b  26 Feb 2019

2.86137e+06 op/sec in 0.37 sec for x^2^2^20
cycles for one square mod N: 771

Estimated time to complete the next test is 11 seconds
2.86068e+06 op/sec in 11.73 sec for x^2^2^25
cycles for one square mod N: 771

x^2^2^38 will take approximately 1 day and 161.474 min (1.11213 day) to compute on this system (t=38)
```

### Command-line tool `sqr`

Here is the sample output of this tool: 
```
./sqr 4 374cd38778e61027a78a9a6f98753f272ca8cbc23b909a8ab041f5030b17abfa
0dcba5a78474cdaff36ae11017d2572497fb4efd1eeeeb560b345c782169cbe09e8bceb30a692ee773f8b9183a87af5d59c112fa685d4edbbf07c5963b5f1ae5bc27548e23b7936d26b432d31cedeef026431499ac7d1a17838c15e03931f956d22ff6896e3e05b60ddbc538ac5173f93d490f4dbc03ccb6470758a34dcafb31285c413d1638417889dcea5c992b46ef4df53d3a3628ed8d62921be794322780f847863b1548e74c6466f9f8614ad68a487f7b10f3df51f0a40c3774b29b7556adf298d7c4585ce2a6b082f1823eb7ddac803d8c12ddfaf3500d8952f288370cfcb28f7753afc8ecc03fe1538fb57e170592b19f2d575c3e933325cbbfa3ea55
```

## API

### sqr\_allocate\_state
Takes three arguments:
- `N` - the modulus
- `Nl` - the size of the    modulus
- `state` - the created object

The function initializes the state, doing all needed pre-calculations that are specific to `N`. 

### sqr_calculate
Takes five arguments:
- `state` - the object, created by `sqr_allocate_state`
- `T` - the number of squaring to perform; usually this is some power of 2, such that `T` = `2^t`
- `xhex` - input value as a string in hexadecial representation
- `outhex` - output buffer
- `outl` - output buffer size in bytes

The function calculates `x^2^T mod N`, returning results as a hexadecimal '\0'-terminated string that is exactly `Nl*2` characters long.

### sqr_free
Takes one arguments:
- `state` - the object to free

The call frees the object `state` that was allocates earlier by `sqr_allocate_state`

## Style

The source code is automatically formatted as follows:
```
make style
```

[uncrustify](https://github.com/uncrustify/uncrustify) tool is needed for the above. 

## License

Standard Apache 2.0 license, matching the license of OpenSSL.


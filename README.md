# Fast Modulo Squaring for  _Verifiable Delay Functions_, _Benchmarking_, and more
> The goal of this project is to provide easy to use fast native code to perform modulo squaring on x86_64.

The project currently offers 4 tools:
1. Squaring API, consisting of 3 APIs in the form Init, Calculate, Free. This API is suitable for integration into other project code base.
2. `sqr_test` tool that provides:
   - Benchmarking of modulo square
   - Report on x86_64 CPU features that are important to achieve best performance
   - Suggests the VDF delay parameter `t`, such that it will take 1 day on the test machine to finish the calculation with this `t`
5. `sqr` command line tool that allows calculation. This is currently the easiest way to interface with this project. 
6. short example on how to integrate `sqr` with JavaScript

2. The reference Prover code written in JavaScript that runs during self-tests.

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

## Building

After cloning the project and changing directory into the location of this [README](README.md):
```
make
```

## Installation

### Dependencies

This project relies on low-level OpenSSL library `libcrypto.so`. On Linux Fedora this library can be installed as follows:
```
dnf install openssl-libs
```
We don't depend on the larger OpenSSL package that includes the `openssl` application and other libraries. 

If you plan to submit code to this respository, you should install `uncrustify` from [uncrustify](https://github.com/uncrustify/uncrustify) somewhere in your `PATH`.

### Componenets

1. The library, consisting of the files `sqrlib.c` and `sqrlib.h`.
2. `sqr` command line tool to perform production-time calculation
3. `sqr_test`, the benchmarking and testing tool, also invoked as 
   ```
   make test
   ```
4. JavaScript [example](example/vdf-prover.js)

## Usage

### Command-line tool `sqr_test`

Here is the sample output of this tool: 
```
./sqr_test
CPU information:
  Brand           : GenuineIntel
  Brand           : Intel(R) Core(TM) i5-3550 CPU @ 3.30GHz
  Features        : ebx=02100800 ecx=7fbae3ff edx=bfebfbff
  Ext. features   : ebx=00000281 ecx=00000000 edx=9c000400
  SSE2            : yes
  AVX             : yes
  MULX, ADX       : no (IMPORTANT!)

OpenSSL version: OpenSSL 1.1.0i-fips  14 Aug 2018

707402 op/sec in 1.48 sec for x^2^2^20
cycles for one square mod N: 4654

Estimated time to complete the next test is 47 seconds
707625 op/sec in 47.42 sec for x^2^2^25
cycles for one square mod N: 4652

x^2^2^36 will take approximately 1 day and 178.547 min (1.12399 day) to compute on this system (t=36)
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

The function initializes the state, doing all needed pre-calculations specific to `N`. 

### sqr_calculate
Takes five arguments:
- `state` - the object, created by sqr\_allocate\_state
- `T` - the number of squaring to perform; usually this is some power of 2, such that `T` = `2^t`
- `xhex` - input value as a string in hexadecial representation
- `outhex` - output buffer
- `outl` - output buffer size in bytes

The function calculates `x^2^T mod N`, returning results as hexadecimal '\0'-tetrminates string that is exactly `Nl*2` characters long.

### sqr_free
Takes one arguments:
- `state` - the object to free

The call frees the object `state` that was allocates earlier by `sqr\_allocate\_state`

## Style

The source code is automatically formatted as follows:
```
make style
```

[uncrustify](https://github.com/uncrustify/uncrustify) tool is needed for the above. 

## License

Standard Apache 2.0 license, matching the license of OpenSSL.

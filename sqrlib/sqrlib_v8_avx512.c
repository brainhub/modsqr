/*
 * Implements v8_square_2048() for 2048 bit modulus via AVX512.
 */
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include <immintrin.h>
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <emmintrin.h>

#include "sqrlib.h"

#include "sqrpriv.h"    // for wbits, W, TRACE, etc

#include "sqrlib_v8_avx512_modN.h"

// use shorter names in this file
#define W v8_W
#define wbits v8_wbits
#define Nl v8_Nl
#define Mw (W - 1)

#define FLOOR8(x) (((x) / 8) * 8)
#define M(x) ((1 << (x)) - 1)

typedef unsigned __int128 uint128_t;

/*
 * Return an array [r,s]: a*b + t = s*W + r
 *
 * This cannot overflow because the max value of a*b + t is (W-1)^2+(W-1) = W^2-W
 */
static inline uint128_t muladdw(uint64_t a, uint64_t b, uint64_t t) {
    return (uint128_t)((uint128_t)a * (uint128_t)b + (uint128_t)t);
}
/*
 * Return 1 if a>b or a==b
 * Return 0 is a<b
 */
static unsigned cmp(const uint64_t *a, const uint64_t *b, unsigned l) {
    l--;
    while(l > 0 && a[l] == b[l])
        l--;
    return a[l] >= b[l];
}
// a -= b
// Assume that each a[i],b[i] are < W
static unsigned sub(uint64_t *a, const uint64_t *b, unsigned l) {
    unsigned i;
    unsigned t;
    unsigned c = 0;
    for( i = 0; i < l; i++) {
        assert(a[i] < W);
        t = a[i] - (b[i] + c);
        c = t > a[i];
        a[i] = t & (W - 1);
    }
    assert(c == 0);
    return c;
}


#if 0
// Return the final carry
static void normalize_test(const uint64_t *a, unsigned l) {
    for(unsigned i = 0; i < l; i++ ) {
        //assert( (a[i] & ~Mw) == 0 );  // actually, this one works too
        assert((a[i] & ~((Mw << 1) | 1)) == 0);
    }
}
#endif


void v8_from_mont(const uint64_t *a, const uint64_t *N, uint64_t Nprime, unsigned l, uint64_t *a_out) {
    uint64_t H, m;
    unsigned i, j;
    uint128_t r;

    memcpy(a_out, a, sizeof(*a) * l);

    for(i = 0; i < l; i++) {
        m = (a_out[0] * Nprime) & (W - 1);
        r = muladdw(m, N[0], a_out[0]);
        H = r >> wbits;

        // ab += N * m
        for(j = 0; j < l - 1; j++) {
            r = muladdw(m, N[j + 1], a_out[j + 1] + H);
            a_out[j] = r & (W - 1);
            H = r >> wbits;
        }

        a_out[l - 1] = H;
    }

    // Last cleanup. Extremely rare, so no need to optimize.
    // I was only able to test this with x = N (which is illegal).
    // Extensive Montecarlo simulations with random x yielded no hit here for wbits=43.
    // We are here for wbits=49 if we pad to v8: we get 0 >= 0.
    if( a_out[l - 1] >= N[l - 1] ) {
        assert(0);
        if(cmp(a_out, N, l) > 0) {
            assert(0);
            sub(a_out, N, l);
        }
    }
}

// a^2 --> out mod N
// Nprime is N': N' * N[0] = -1 mod W
// All buffers are 64-byte aligned.
// Sizes:
//   a, N : Nl (48) elements, zero padded if needed
//   Nprime: wbits
//   out: 7+2*Nl, where top 7 and bottom v8_Wl elements are zeroes.
//   The out has this condition on entrance as well (and it's maintained on exit).

void v8_square_2048(uint64_t T, const uint64_t a[v8_Nl], const uint64_t N[v8_Nl], uint64_t Nprime, uint64_t out[7 + 2 * v8_Nl]) {
    const __m512i zero = _mm512_setzero_si512();

    __m512i O_5 =  _mm512_load_si512(a + 0);
    __m512i O_6 = _mm512_load_si512(a + 1 * 8);
    __m512i O_7 = _mm512_load_si512(a + 2 * 8);
    __m512i O_8 = _mm512_load_si512(a + 3 * 8);
    __m512i O_9 = _mm512_load_si512(a + 4 * 8);
    __m512i O_10 = _mm512_load_si512(a + 5 * 8);

    while(T--) {

        __m512i O_2;
        __m512i O_3;
        __m512i O_4;

        __m512i N8_0 = _mm512_load_si512(N + 0);
        __m512i N8_1 = _mm512_load_si512(N + 1 * 8);
        __m512i N8_2 = _mm512_load_si512(N + 2 * 8);
        __m512i N8_3 = _mm512_load_si512(N + 3 * 8);
        __m512i N8_4 = _mm512_load_si512(N + 4 * 8);
        __m512i N8_5 = _mm512_load_si512(N + 5 * 8);

        // i below is an iteration of the outer loop, i=0...40

        { // loop over A
            const __m512i A8_1 = O_6;
            const __m512i A8_2 = O_7;
            const __m512i A8_3 = O_8;
            const __m512i A8_4 = O_9;
            const __m512i A8_5 = O_10;

            { // i = 0,1,2,3,4,5,6
                const __m512i A8_0 = O_5;
                __m512i Aj8_0, Aj8_1, Aj8_2, Aj8_3, Aj8_4, Aj8_5;
                __m512i Ai8;

                const __m512i A8_0_sqr  = _mm512_slli_epi64(A8_0, 1);// square once, then reuse

                __m512i O_0;
                __m512i Oh_0;
                __m512i O_1;
                __m512i Oh_1;

                {
                    const __m512i _0123 = _mm512_set_epi64(3,3, 2,2, 1,1, 0,0);
                    const __m512i _4567 = _mm512_set_epi64(7,7, 6,6, 5,5, 4,4);
                    __m512i A8_0_0 = _mm512_permutexvar_epi64(_0123, A8_0);
                    __m512i A8_0_1 = _mm512_permutexvar_epi64(_4567, A8_0);
                    __m512i t2 = _mm512_permutexvar_epi64(_0123, A8_5);
                    O_0 = _mm512_mask_madd52lo_epu64(zero, 0x55, A8_0_0, A8_0_0); // start
                    Oh_0 = _mm512_mask_madd52hi_epu64(zero, 0xaa, A8_0_0, A8_0_0); //start
                    O_1 = _mm512_mask_madd52lo_epu64(zero, 0x55, A8_0_1, A8_0_1); // start
                    Oh_1 = _mm512_mask_madd52hi_epu64(zero, 0xaa, A8_0_1, A8_0_1); // start
                    O_10 = _mm512_mask_madd52lo_epu64(zero, 3 & 0x55, t2, t2);        // 1 elements out of 2
                }

                __m512i Oh_2;
                __m512i Oh_3;
                __m512i Oh_4;
                __m512i Oh_5;
                __m512i Oh_6;

                // i==0
                Ai8 = _mm512_permutexvar_epi64(zero, A8_0_sqr); // [A0[0],A0[0], ... A0[0]]
                Aj8_0 = _mm512_maskz_alignr_epi64(M(7) << 1, A8_0, A8_0,  8 - 1);  // sAj8 = A8 << 64
                Aj8_1 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 1);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 1);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 1);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 1);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 1);

                //LOAD_6(0*8+0);

                //MUL_6();

                O_0 = _mm512_mask_madd52lo_epu64(O_0, 0xff ^ M(1), Ai8, A8_0); // element A8_0[0] is dropped, shift zero
                O_1 = _mm512_madd52lo_epu64(O_1, Ai8, A8_1);
                O_2 = _mm512_madd52lo_epu64(zero, Ai8, A8_2); // start
                O_3 = _mm512_madd52lo_epu64(zero, Ai8, A8_3); // start
                O_4 = _mm512_madd52lo_epu64(zero, Ai8, A8_4); // start
                O_5 = _mm512_madd52lo_epu64(zero, Ai8, A8_5); // start

                Oh_0 = _mm512_mask_madd52hi_epu64(Oh_0, 0xff ^ M(2), Ai8, Aj8_0); // two elements less, shift one
                Oh_1 = _mm512_madd52hi_epu64(Oh_1, Ai8, Aj8_1);
                Oh_2 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_2); // start
                Oh_3 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_3); // start
                Oh_4 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_4); // start
                Oh_5 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_5); // start

                //STORE_6(0*8+0);

                // i == 1
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(1), A8_0_sqr); // [A0[1],A0[1], ... A0[1]]
                //LOAD_6(0*8+1);
                //MUL_6();
                O_0 = _mm512_mask_madd52lo_epu64(O_0, 0xff ^ M(3), Ai8, Aj8_0); // 3 elelents less, shift one
                O_1 = _mm512_madd52lo_epu64(O_1, Ai8, Aj8_1);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_2);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_3);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_4);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_5);

                Aj8_0 = _mm512_maskz_alignr_epi64(M(8 - 4) << 4, A8_0, A8_0,  8 - 2);// << 64*2
                Aj8_1 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 2);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 2);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 2);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 2);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 2);

                Oh_0 = _mm512_madd52hi_epu64(Oh_0, Ai8, Aj8_0);
                Oh_1 = _mm512_madd52hi_epu64(Oh_1, Ai8, Aj8_1);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_2);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_3);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_4);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_5);


                // i==2
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(2), A8_0_sqr);
                //LOAD_6(0*8+2);
                //MUL_6();
                O_0 = _mm512_mask_madd52lo_epu64(O_0, 0xff ^ M(5), Ai8, Aj8_0);
                O_1 = _mm512_madd52lo_epu64(O_1, Ai8, Aj8_1);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_2);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_3);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_4);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_5);

                Aj8_0 = _mm512_maskz_alignr_epi64(M(8 - 6) << 6, A8_0, A8_0,  8 - 3);// << 64*3
                Aj8_1 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 3);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 3);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 3);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 3);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 3);

                Oh_0 = _mm512_madd52hi_epu64(Oh_0, Ai8, Aj8_0); // last
                Oh_1 = _mm512_madd52hi_epu64(Oh_1, Ai8, Aj8_1);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_2);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_3);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_4);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_5);
                //STORE_6(0*8+2);

                // i==3 (jump the Oh)
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(3), A8_0_sqr);
                //LOAD_6(0*8+3);
                //MUL_6();
                O_0 = _mm512_mask_madd52lo_epu64(_mm512_add_epi64(O_0, _mm512_slli_epi64(Oh_0, (52 - wbits))),  0xff ^ M(7),Ai8, Aj8_0);// last
                O_1 = _mm512_madd52lo_epu64(O_1, Ai8, Aj8_1);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_2);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_3);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_4);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_5);

                Aj8_1 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 4);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 4);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 4);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 4);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 4);

                Oh_1 = _mm512_madd52hi_epu64(Oh_1, Ai8, Aj8_1);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_2);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_3);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_4);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_5);
                //STORE_6(0*8+3);

                // i== 4
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(4), A8_0_sqr);
                //LOAD_6(0*8+4);

                //MUL_6();
                O_1 = _mm512_mask_madd52lo_epu64(O_1, 0xff ^ M(1), Ai8, Aj8_1);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_2);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_3);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_4);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_5);

                Aj8_0 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 5);
                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 5);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 5);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 5);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 5);

                Oh_1 = _mm512_mask_madd52hi_epu64(Oh_1, 0xff ^ M(2), Ai8, Aj8_0);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);

                // i==5
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(5), A8_0_sqr);

                O_1 = _mm512_mask_madd52lo_epu64(O_1, 0xff ^ M(3), Ai8, Aj8_0);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);

                Aj8_0 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 6);
                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 6);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 6);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 6);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 6);

                Oh_1 = _mm512_mask_madd52hi_epu64(Oh_1,  0xff ^ M(4), Ai8, Aj8_0);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);

                // i==6
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(6), A8_0_sqr);

                O_1 = _mm512_mask_madd52lo_epu64(O_1, 0xff ^ M(5), Ai8, Aj8_0);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);

                Aj8_0 = _mm512_alignr_epi64(A8_1, A8_0,  8 - 7);
                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 7);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 7);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 7);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 7);

                Oh_1 = _mm512_mask_madd52hi_epu64(Oh_1,  0xff ^ M(6), Ai8, Aj8_0);
                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);

                //assert(i == 7);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(7), A8_0_sqr);

                O_1 = _mm512_mask_madd52lo_epu64(O_1, 0xff ^ M(7), Ai8, Aj8_0);
                O_2 = _mm512_madd52lo_epu64(O_2, Ai8, Aj8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);

                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, A8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, A8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, A8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, A8_4);
                Oh_6 = _mm512_madd52hi_epu64(zero, Ai8, A8_5);

                Oh_0 = zero;
                O_7 = zero;

                // Mid-way normalize. Needed for wbits=49 and higher, but not lower
                {
                    const __m512i WMASK = _mm512_set1_epi64(W - 1);
                    __m512i O_2c, O_3c, O_4c, O_5c;

                    O_2c = _mm512_srli_epi64(O_2, wbits);
                    O_3c = _mm512_srli_epi64(O_3, wbits);
                    O_4c = _mm512_srli_epi64(O_4, wbits);
                    O_5c = _mm512_srli_epi64(O_5, wbits);

                    O_2 = _mm512_add_epi64(_mm512_and_epi64(O_2, WMASK), _mm512_alignr_epi64(O_2c, zero, 8 - 1));
                    O_3 = _mm512_add_epi64(_mm512_and_epi64(O_3, WMASK), _mm512_alignr_epi64(O_3c, O_2c, 8 - 1));
                    O_4 = _mm512_add_epi64(_mm512_and_epi64(O_4, WMASK), _mm512_alignr_epi64(O_4c, O_3c, 8 - 1));
                    O_5 = _mm512_add_epi64(_mm512_and_epi64(O_5, WMASK), _mm512_alignr_epi64(O_5c, O_4c, 8 - 1));
                    O_6 = _mm512_maskz_alignr_epi64(M(1), O_5c, O_5c, 8 - 1);
                }

                RED_N(0, O_0,O_1,O_2,O_3,O_4,O_5,O_6, Oh_0,Oh_1,Oh_2,Oh_3,Oh_4,Oh_5);
                RED_N(8, O_1,O_2,O_3,O_4,O_5,O_6,O_7, Oh_1,Oh_2,Oh_3,Oh_4,Oh_5,Oh_6);

                // Add carries now to save registers. This is a few cycles faster than doing this at the end
                O_2 = _mm512_add_epi64(O_2, _mm512_slli_epi64(Oh_2, 52 - wbits));
                O_3 = _mm512_add_epi64(O_3, _mm512_slli_epi64(Oh_3, 52 - wbits));
                O_4 = _mm512_add_epi64(O_4, _mm512_slli_epi64(Oh_4, 52 - wbits));
                O_5 = _mm512_add_epi64(O_5, _mm512_slli_epi64(Oh_5, 52 - wbits));
                O_6 = _mm512_add_epi64(O_6, _mm512_slli_epi64(Oh_6, 52 - wbits));
            }

            { // i== 7,8,9,10,11,12,13,14
                __m512i Aj8_1, Aj8_2, Aj8_3, Aj8_4, Aj8_5;
                __m512i Ai8;
                const __m512i A8_1_sqr  = _mm512_slli_epi64(A8_1, 1);// square once, then reuse

                // these Oh_i are assumed zero after above
                __m512i Oh_2;
                __m512i Oh_3;
                __m512i Oh_4;
                __m512i Oh_5;
                __m512i Oh_6;
                __m512i Oh_7;
                __m512i Oh_8;

                const __m512i _0123 = _mm512_set_epi64(3,3, 2,2, 1,1, 0,0);
                const __m512i _4567 = _mm512_set_epi64(7,7, 6,6, 5,5, 4,4);
                {
                    __m512i t = _mm512_permutexvar_epi64(_0123, A8_1);
                    O_2 = _mm512_mask_madd52lo_epu64(O_2, 0x55, t, t);
                    Oh_2 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }
                {
                    __m512i t = _mm512_permutexvar_epi64(_4567, A8_1);
                    O_3 = _mm512_mask_madd52lo_epu64(O_3, 0x55, t, t);
                    Oh_3 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }
                {
                    __m512i t = _mm512_permutexvar_epi64(_0123, A8_2);
                    O_4 = _mm512_mask_madd52lo_epu64(O_4, 0x55, t, t);
                    Oh_4 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }
                {
                    __m512i t = _mm512_permutexvar_epi64(_4567, A8_2);
                    O_5 = _mm512_mask_madd52lo_epu64(O_5, 0x55, t, t);
                    Oh_5 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }

                // i == 8
                Ai8 = _mm512_permutexvar_epi64(zero, A8_1_sqr); // [A1[0],A1[0], ... A1[0]]

                O_2 = _mm512_mask_madd52lo_epu64(O_2, 0xff ^ M(1), Ai8, A8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, A8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, A8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, A8_4);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, A8_5);

                Aj8_1 = _mm512_maskz_alignr_epi64(M(8 - 2) << 2, A8_1, A8_1, 8 - 1);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 1);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 1);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 1);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 1);

                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);
                Oh_6 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_5);

                // assert(i == 8+1);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(1), A8_1_sqr); // [A1[1],A1[1], ... A1[1]]

                O_2 = _mm512_mask_madd52lo_epu64(O_2, 0xff ^ M(3), Ai8, Aj8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_5);

                Aj8_1 = _mm512_maskz_alignr_epi64(M(8 - 4) << 4, A8_1, A8_1, 8 - 2);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 2);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 2);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 2);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 2);

                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1);
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_5);

                //assert(i == 8+2);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(2), A8_1_sqr);
                O_2 = _mm512_mask_madd52lo_epu64(O_2, 0xff ^ M(5), Ai8, Aj8_1);
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_5);

                Aj8_1 = _mm512_maskz_alignr_epi64(M(8 - 6) << 6, A8_1, A8_1, 8 - 3);
                Aj8_2 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 3);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 3);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 3);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 3);

                Oh_2 = _mm512_madd52hi_epu64(Oh_2, Ai8, Aj8_1); // last
                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_2);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_3);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_4);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_5);

                //assert(i == 8+3);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(3), A8_1_sqr);
                O_2 = _mm512_mask_madd52lo_epu64(O_2, 0xff ^ M(7), Ai8, Aj8_1); // last
                O_3 = _mm512_madd52lo_epu64(O_3, Ai8, Aj8_2);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_3);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_4);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_5);

                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 4);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 4);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 4);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 4);

                Oh_3 = _mm512_madd52hi_epu64(Oh_3, Ai8, Aj8_1);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);

                //assert(i == 8+4);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(4), A8_1_sqr);
                O_3 = _mm512_mask_madd52lo_epu64(O_3, 0xff ^ M(1), Ai8, Aj8_1);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);

                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 5);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 5);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 5);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 5);

                Oh_3 = _mm512_mask_madd52hi_epu64(Oh_3, 0xff ^ M(2), Ai8, Aj8_1);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);

                //assert(i == 8+5);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(5), A8_1_sqr);
                O_3 = _mm512_mask_madd52lo_epu64(O_3, 0xff ^ M(3), Ai8, Aj8_1);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);

                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 6);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 6);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 6);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 6);

                Oh_3 = _mm512_mask_madd52hi_epu64(Oh_3, 0xff ^ M(4), Ai8, Aj8_1);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);

                //assert(i == 8+6);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(6), A8_1_sqr);

                O_3 = _mm512_mask_madd52lo_epu64(O_3, 0xff ^ M(5), Ai8, Aj8_1);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);

                Aj8_1 = _mm512_alignr_epi64(A8_2, A8_1,  8 - 7);
                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 7);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 7);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 7);

                Oh_3 = _mm512_mask_madd52hi_epu64(Oh_3, 0xff ^ M(6), Ai8, Aj8_1);
                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);

                //assert(i == 1*8 + 7);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(7), A8_1_sqr);

                O_3 = _mm512_mask_madd52lo_epu64(O_3, 0xff ^ M(7), Ai8, Aj8_1);
                O_4 = _mm512_madd52lo_epu64(O_4, Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);

                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, A8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, A8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, A8_4);
                Oh_7 = _mm512_madd52hi_epu64(zero, Ai8, A8_5);

                O_8 = zero;
                O_9 = zero;
                Oh_8 = zero;
                RED_N(2 * 8, O_2,O_3,O_4,O_5,O_6,O_7,O_8, Oh_2,Oh_3,Oh_4,Oh_5,Oh_6,Oh_7);
                RED_N(3 * 8, O_3,O_4,O_5,O_6,O_7,O_8,O_9, Oh_3,Oh_4,Oh_5,Oh_6,Oh_7,Oh_8);

                O_4 = _mm512_add_epi64(O_4, _mm512_slli_epi64(Oh_4, 52 - wbits));
                O_5 = _mm512_add_epi64(O_5, _mm512_slli_epi64(Oh_5, 52 - wbits));
                O_6 = _mm512_add_epi64(O_6, _mm512_slli_epi64(Oh_6, 52 - wbits));
                O_7 = _mm512_add_epi64(O_7, _mm512_slli_epi64(Oh_7, 52 - wbits));
                O_8 = _mm512_add_epi64(O_8, _mm512_slli_epi64(Oh_8, 52 - wbits));
                // done with O_3, Oh_3
            }

            { // i = 15,16,17,18,19,20,21,22
                __m512i Aj8_2, Aj8_3, Aj8_4, Aj8_5;
                const __m512i A8_2_sqr  = _mm512_slli_epi64(A8_2, 1);// square once, then reuse
                __m512i Ai8;

                __m512i Oh_4;
                __m512i Oh_5;
                __m512i Oh_6;
                __m512i Oh_7;
                __m512i Oh_8;
                __m512i Oh_9;
                __m512i Oh_10;

                {
                    const __m512i _0123 = _mm512_set_epi64(3,3, 2,2, 1,1, 0,0);
                    //const __m512i _4567 = _mm512_set_epi64(7,7, 6,6, 5,5, 4,4);
                    //__m512i t = _mm512_permutexvar_epi64(_4567, A8_4);
                    __m512i t2 = _mm512_permutexvar_epi64(_0123, A8_5);

                    Oh_10 = _mm512_mask_madd52hi_epu64(zero, 3 & 0xaa, t2, t2);
                }

                //assert(i == 2*8 + 0 );
                Ai8 = _mm512_permutexvar_epi64(zero, A8_2_sqr); // [A2[0],A2[0], ... A2[0]]

                Aj8_2 = _mm512_maskz_alignr_epi64(M(8 - 2) << 2, A8_2, A8_2, 8 - 1);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2, 8 - 1);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3, 8 - 1);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4, 8 - 1);

                O_4 = _mm512_mask_madd52lo_epu64(O_4, 0xff ^ M(1), Ai8, A8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, A8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, A8_4);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, A8_5);

                Oh_4 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_4);
                Oh_7 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_5);

                //assert(i == 2*8 + 1);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(1), A8_2_sqr); // [A2[1],A2[1], ... A2[1]]
                O_4 = _mm512_mask_madd52lo_epu64(O_4, 0xff ^ M(3), Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_5);

                Aj8_2 = _mm512_maskz_alignr_epi64(M(8 - 4) << 4, A8_2, A8_2, 8 - 2);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2, 8 - 2);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3, 8 - 2);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4, 8 - 2);

                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_5);

                //assert(i == 2*8 + 2);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(2), A8_2_sqr);
                O_4 = _mm512_mask_madd52lo_epu64(O_4, 0xff ^ M(5), Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_5);

                Aj8_2 = _mm512_maskz_alignr_epi64(M(8 - 6) << 6, A8_2, A8_2, 8 - 3);
                Aj8_3 = _mm512_alignr_epi64(A8_3, A8_2, 8 - 3);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3, 8 - 3);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4, 8 - 3);

                Oh_4 = _mm512_madd52hi_epu64(Oh_4, Ai8, Aj8_2);
                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_3);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_4);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_5);

                //assert(i == 2*8 + 3);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(3), A8_2_sqr);
                O_4 = _mm512_mask_madd52lo_epu64(O_4, 0xff ^ M(7), Ai8, Aj8_2);
                O_5 = _mm512_madd52lo_epu64(O_5, Ai8, Aj8_3);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_4);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_5);

                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 4);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 4);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 4);

                Oh_5 = _mm512_madd52hi_epu64(Oh_5, Ai8, Aj8_2);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);

                //assert(i == 2*8 + 4);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(4), A8_2_sqr);
                O_5 = _mm512_mask_madd52lo_epu64(O_5, 0xff ^ M(1), Ai8, Aj8_2);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);

                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 5);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 5);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 5);

                Oh_5 = _mm512_mask_madd52hi_epu64(Oh_5, 0xff ^ M(2), Ai8, Aj8_2);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);

                //assert(i == 2*8 + 5);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(5), A8_2_sqr);
                O_5 = _mm512_mask_madd52lo_epu64(O_5, 0xff ^ M(3), Ai8, Aj8_2);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);

                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 6);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 6);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 6);

                Oh_5 = _mm512_mask_madd52hi_epu64(Oh_5, 0xff ^ M(4), Ai8, Aj8_2);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);

                //assert(i == 2*8 + 6);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(6), A8_2_sqr);
                O_5 = _mm512_mask_madd52lo_epu64(O_5, 0xff ^ M(5), Ai8, Aj8_2);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);

                Aj8_2 = _mm512_alignr_epi64(A8_3, A8_2,  8 - 7);
                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 7);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 7);

                Oh_5 = _mm512_mask_madd52hi_epu64(Oh_5, 0xff ^ M(6), Ai8, Aj8_2);
                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);

                //assert(i == 2*8+7);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(7), A8_2_sqr);
                O_5 = _mm512_mask_madd52lo_epu64(O_5, 0xff ^ M(7), Ai8, Aj8_2);
                O_6 = _mm512_madd52lo_epu64(O_6, Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);

                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, A8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, A8_4);
                Oh_8 = _mm512_madd52hi_epu64(zero, Ai8, A8_5);

                {
                    const __m512i _4567 = _mm512_set_epi64(7,7, 6,6, 5,5, 4,4);
                    __m512i t = _mm512_permutexvar_epi64(_4567, A8_4);
                    Oh_9 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                    O_9 = _mm512_mask_madd52lo_epu64(O_9, 0x55, t, t);
                }

                RED_N(4 * 8, O_4,O_5,O_6,O_7,O_8,O_9,O_10, Oh_4,Oh_5,Oh_6,Oh_7,Oh_8,Oh_9);

                // finally, reduce by one element, for v8_Nl_tight in total
                RED_N_1(5 * 8, O_5,O_6,O_7,O_8,O_9,O_10, Oh_5,Oh_6,Oh_7,Oh_8,Oh_9,Oh_10);

                // done with O_4, Oh_4

                O_5 = _mm512_add_epi64(O_5, _mm512_slli_epi64(Oh_5, 52 - wbits));
                O_6 = _mm512_add_epi64(O_6, _mm512_slli_epi64(Oh_6, 52 - wbits));
                O_7 = _mm512_add_epi64(O_7, _mm512_slli_epi64(Oh_7, 52 - wbits));
                O_8 = _mm512_add_epi64(O_8, _mm512_slli_epi64(Oh_8, 52 - wbits));
                O_9 = _mm512_add_epi64(O_9, _mm512_slli_epi64(Oh_9, 52 - wbits));
                O_10 = _mm512_add_epi64(O_10, _mm512_slli_epi64(Oh_10, 52 - wbits));
            }

            // i = 23,24,25,26, 27,28,29,30
            {
                __m512i Aj8_3, Aj8_4, Aj8_5;
                const __m512i A8_3_sqr  = _mm512_slli_epi64(A8_3, 1);// square once, then reuse
                __m512i Ai8;

                __m512i Oh_6;
                __m512i Oh_7;
                __m512i Oh_8;
                __m512i Oh_9;

                const __m512i _0123 = _mm512_set_epi64(3,3, 2,2, 1,1, 0,0);
                const __m512i _4567 = _mm512_set_epi64(7,7, 6,6, 5,5, 4,4);
                {
                    __m512i t = _mm512_permutexvar_epi64(_0123, A8_3);
                    O_6 = _mm512_mask_madd52lo_epu64(O_6, 0x55, t, t);
                    Oh_6 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }
                {
                    __m512i t = _mm512_permutexvar_epi64(_4567, A8_3);
                    O_7 = _mm512_mask_madd52lo_epu64(O_7, 0x55, t, t);
                    Oh_7 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }
                {
                    __m512i t = _mm512_permutexvar_epi64(_0123, A8_4);
                    O_8 = _mm512_mask_madd52lo_epu64(O_8, 0x55, t, t);
                    Oh_8 = _mm512_mask_madd52hi_epu64(zero, 0xaa, t, t);
                }

                //assert(i == 3*8);
                Ai8 = _mm512_permutexvar_epi64(zero, A8_3_sqr); // [A3[0],A3[0], ... A3[0]]
                O_6 = _mm512_mask_madd52lo_epu64(O_6, 0xff ^ M(1), Ai8, A8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, A8_4);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, A8_5);

                Aj8_3 = _mm512_maskz_alignr_epi64(M(8 - 2) << 2, A8_3, A8_3, 8 - 1);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 1);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 1);

                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_5);

                //assert(i == 3*8+1);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(1), A8_3_sqr); // [A3[1],A3[1], ... A3[1]]
                O_6 = _mm512_mask_madd52lo_epu64(O_6, 0xff ^ M(3), Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_5);

                Aj8_3 = _mm512_maskz_alignr_epi64(M(8 - 4) << 4, A8_3, A8_3, 8 - 2);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 2);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 2);

                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_5);

                //assert(i == 3*8+2);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(2), A8_3_sqr);
                O_6 = _mm512_mask_madd52lo_epu64(O_6, 0xff ^ M(5), Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_5);

                Aj8_3 = _mm512_maskz_alignr_epi64(M(8 - 6) << 6, A8_3, A8_3, 8 - 3);
                Aj8_4 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 3);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 3);

                Oh_6 = _mm512_madd52hi_epu64(Oh_6, Ai8, Aj8_3);
                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_4);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_5);

                //assert(i == 3*8+3);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(3), A8_3_sqr);
                O_6 = _mm512_mask_madd52lo_epu64(O_6, 0xff ^ M(7), Ai8, Aj8_3);
                O_7 = _mm512_madd52lo_epu64(O_7, Ai8, Aj8_4);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_5);

                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 4);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 4);

                Oh_7 = _mm512_madd52hi_epu64(Oh_7, Ai8, Aj8_3);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);

                //assert(i == 3*8+4);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(4), A8_3_sqr);
                O_7 = _mm512_mask_madd52lo_epu64(O_7, 0xff ^ M(1), Ai8, Aj8_3);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_4);

                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 5);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 5);

                Oh_7 = _mm512_mask_madd52hi_epu64(Oh_7, 0xff ^ M(2), Ai8, Aj8_3);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);

                //assert(i == 3*8+5);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(5), A8_3_sqr);
                O_7 = _mm512_mask_madd52lo_epu64(O_7, 0xff ^ M(3), Ai8, Aj8_3);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_4);

                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 6);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 6);

                Oh_7 = _mm512_mask_madd52hi_epu64(Oh_7, 0xff ^ M(4), Ai8, Aj8_3);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);

                //assert(i == 3*8+6);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(6), A8_3_sqr);
                O_7 = _mm512_mask_madd52lo_epu64(O_7, 0xff ^ M(5), Ai8, Aj8_3);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_4);

                Aj8_3 = _mm512_alignr_epi64(A8_4, A8_3,  8 - 7);
                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 7);

                Oh_7 = _mm512_mask_madd52hi_epu64(Oh_7, 0xff ^ M(6), Ai8, Aj8_3);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);

                //assert(i == 3*8+7);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(7), A8_3_sqr);
                O_7 = _mm512_mask_madd52lo_epu64(O_7, 0xff ^ M(7), Ai8, Aj8_3);
                O_8 = _mm512_madd52lo_epu64(O_8, Ai8, Aj8_4);
                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, A8_4);
                Oh_9 = _mm512_madd52hi_epu64(zero, Ai8, A8_5);

                O_6 = _mm512_add_epi64(O_6, _mm512_slli_epi64(Oh_6, 52 - wbits)); // last use of Oh_6
                O_7 = _mm512_add_epi64(O_7, _mm512_slli_epi64(Oh_7, 52 - wbits)); // last use of Oh_7
                O_8 = _mm512_add_epi64(O_8, _mm512_slli_epi64(Oh_8, 52 - wbits));
                O_9 = _mm512_add_epi64(O_9, _mm512_slli_epi64(Oh_9, 52 - wbits));
            }

            { // i = 31,32,33,34,35,36,37,38
                __m512i Aj8_4, Aj8_5, Aj8_6;
                const __m512i A8_4_sqr  = _mm512_slli_epi64(A8_4, 1);// square once, then reuse
                __m512i Ai8;
                __m512i Oh_8;
                __m512i Oh_9;
                __m512i Oh_10;

                //assert(i == 4*8);
                Ai8 = _mm512_permutexvar_epi64(zero, A8_4_sqr); // [A4[0],A4[0], ... A4[0]]
                O_8 = _mm512_mask_madd52lo_epu64(O_8, 0xff ^ M(1), Ai8, A8_4);
                O_9 = _mm512_madd52lo_epu64(O_9, Ai8, A8_5);

                Aj8_4 = _mm512_maskz_alignr_epi64(M(8 - 2) << 2, A8_4, A8_4, 8 - 1);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 1);

                Oh_8 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_4);
                Oh_9 = _mm512_madd52hi_epu64(zero, Ai8, Aj8_5);

                //assert(i == 4*8+1);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(1), A8_4_sqr); // [A4[1],A4[1], ... A4[1]]
                O_8 = _mm512_mask_madd52lo_epu64(O_8, 0xff ^ M(3), Ai8, Aj8_4);
                O_9 = _mm512_madd52lo_epu64(O_9, Ai8, Aj8_5);

                Aj8_4 = _mm512_maskz_alignr_epi64(M(8 - 4) << 4, A8_4, A8_4, 8 - 2);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 2);

                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);
                Oh_9 = _mm512_madd52hi_epu64(Oh_9, Ai8, Aj8_5);

                //assert(i == 4*8+2);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(2), A8_4_sqr); // [A4[1],A4[1], ... A4[1]]
                O_8 = _mm512_mask_madd52lo_epu64(O_8, 0xff ^ M(5), Ai8, Aj8_4);
                O_9 = _mm512_madd52lo_epu64(O_9, Ai8, Aj8_5);

                Aj8_4 = _mm512_maskz_alignr_epi64(M(8 - 6) << 6, A8_4, A8_4, 8 - 3);
                Aj8_5 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 3);

                Oh_8 = _mm512_madd52hi_epu64(Oh_8, Ai8, Aj8_4);
                Oh_9 = _mm512_madd52hi_epu64(Oh_9, Ai8, Aj8_5);

                //assert(i == 4*8+3);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(3), A8_4_sqr);
                O_8 = _mm512_mask_madd52lo_epu64(O_8, 0xff ^ M(7), Ai8, Aj8_4);
                O_9 = _mm512_madd52lo_epu64(O_9, Ai8, Aj8_5);

                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 4);

                Oh_9 = _mm512_madd52hi_epu64(Oh_9, Ai8, Aj8_4);

                //assert(i == 4*8+4);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(4), A8_4_sqr);
                O_9 = _mm512_mask_madd52lo_epu64(O_9, 0xff ^ M(1), Ai8, Aj8_4);

                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 5);

                Oh_9 =  _mm512_mask_madd52hi_epu64(Oh_9, 0xff ^ M(2), Ai8, Aj8_4);

                //assert(i == 4*8+5);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(5), A8_4_sqr);
                O_9 = _mm512_mask_madd52lo_epu64(O_9, 0xff ^ M(3), Ai8, Aj8_4);

                Aj8_4 = _mm512_alignr_epi64(A8_5, A8_4,  8 - 6);

                Oh_9 =  _mm512_mask_madd52hi_epu64(Oh_9, 0xff ^ M(4), Ai8, Aj8_4);


                //assert(i == 4*8+6);
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(6), A8_4_sqr);

                O_9 = _mm512_mask_madd52lo_epu64(O_9, 0xff ^ M(5), Ai8, Aj8_4);

                Oh_9 =  _mm512_mask_madd52hi_epu64(Oh_9, 0xff ^ M(6), Ai8, _mm512_alignr_epi64(A8_5, A8_4,  8 - 7));

                //assert(i == 4*8+7);
                // 2 cycles faster to do this in a naive way on Cannon Lake
                Ai8 = _mm512_permutexvar_epi64(_mm512_set1_epi64(7), A8_4_sqr);
                Aj8_6 = _mm512_alignr_epi64(O_10, O_9, 7);
                Oh_10 = _mm512_madd52hi_epu64(zero, Ai8, A8_5);         // 1 element
                Aj8_6 = _mm512_madd52lo_epu64(Aj8_6, Ai8, A8_5);        // 1 element
                O_8 = _mm512_add_epi64(O_8, _mm512_slli_epi64(Oh_8, 52 - wbits)); // last use of Oh_8
                O_9 =  _mm512_mask_alignr_epi64(O_9, 0xff ^ M(7), Aj8_6, Aj8_6, 8 - (7));

                O_10 = _mm512_add_epi64(O_10, _mm512_slli_epi64(Oh_10, 52 - wbits));
                O_9 = _mm512_add_epi64(O_9, _mm512_slli_epi64(Oh_9, 52 - wbits));       // done with Oh_9

                // i = 40 (5*8) -- nothing
                //assert(i == 5*8);
                //assert(i == 5*8+1);
            }
        } // loop over A

        /* at this point we have elements in O_5, O_6, O_7, O_8, O_9, O_10. O_5 has 7 top, and O_10 has 2 lower elements, for 41 in total */
        O_5 = _mm512_alignr_epi64(O_6, O_5, 1);
        O_6 = _mm512_alignr_epi64(O_7, O_6, 1);
        O_7 = _mm512_alignr_epi64(O_8, O_7, 1);
        O_8 = _mm512_alignr_epi64(O_9, O_8, 1);
        O_9 = _mm512_alignr_epi64(O_10, O_9, 1);
        O_10 = _mm512_maskz_alignr_epi64(1, O_10, O_10, 1);      // O_10 gets 1 element

        // normalize
        {
            const __m512i WMASK = _mm512_set1_epi64(W - 1);
            __m512i O_5c, O_6c, O_7c, O_8c, O_9c, O_10c;

            O_5c = _mm512_srli_epi64(O_5, wbits); // carries of O_5
            O_6c = _mm512_srli_epi64(O_6, wbits);
            O_7c = _mm512_srli_epi64(O_7, wbits);
            O_8c = _mm512_srli_epi64(O_8, wbits);
            O_9c = _mm512_srli_epi64(O_9, wbits);
            O_10c = _mm512_srli_epi64(O_10, wbits);

            O_5 = _mm512_add_epi64(_mm512_and_epi64(O_5, WMASK), _mm512_maskz_alignr_epi64(M(8 - 1) << 1, O_5c, O_5c, 8 - 1)); // O_5 += carries_of_O_5 << 64
            O_6 = _mm512_add_epi64(_mm512_and_epi64(O_6, WMASK), _mm512_alignr_epi64(O_6c, O_5c, 8 - 1)); // O_6 += carries_of_O_6 || carries_of_O5
            O_7 = _mm512_add_epi64(_mm512_and_epi64(O_7, WMASK), _mm512_alignr_epi64(O_7c, O_6c, 8 - 1));
            O_8 = _mm512_add_epi64(_mm512_and_epi64(O_8, WMASK), _mm512_alignr_epi64(O_8c, O_7c, 8 - 1));
            O_9 = _mm512_add_epi64(_mm512_and_epi64(O_9, WMASK), _mm512_alignr_epi64(O_9c, O_8c, 8 - 1));
            O_10 = _mm512_add_epi64(_mm512_and_epi64(O_10, WMASK), _mm512_alignr_epi64(O_10c, O_9c, 8 - 1));
        }

    }   // loop for T

    // final store
    _mm512_store_si512(out + 6 * 8, O_5);
    _mm512_store_si512(out + 7 * 8, O_6);
    _mm512_store_si512(out + 8 * 8, O_7);
    _mm512_store_si512(out + 9 * 8, O_8);
    _mm512_store_si512(out + 10 * 8, O_9);
    _mm512_store_si512(out + 11 * 8, O_10);

    // We could get rid of this clearing, but this doesn't matter for performance.
    _mm512_store_si512(out + 0 * 8, zero);
    _mm512_store_si512(out + 1 * 8, zero);
    _mm512_store_si512(out + 2 * 8, zero);
    _mm512_store_si512(out + 3 * 8, zero);
    _mm512_store_si512(out + 4 * 8, zero);
    _mm512_store_si512(out + 5 * 8, zero);

    //SQR_TRACE_BUF(out, 2*Nl, "output" );
    //normalize_test(out, Nl * 2);
}



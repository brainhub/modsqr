/*
 * Implements sqr_allocate_state() / sqr_calculate() / sqr_free_state()
 *
 * Depending on CPU capability, this takes OpenSSL code path or OpenSSL+AVX512 optimized code.
 *
 */
#include <stdio.h>
#include <openssl/bn.h>
#include <assert.h>
#include <string.h>
#include <cpuid.h>

#include "sqrlib.h"
#include "sqrpriv.h"


typedef struct Sqr_State {
    BIGNUM *N;
    unsigned Nl;
    BN_MONT_CTX *montN;         // Montgomery multiplication context for mod N
    BN_CTX *bnctx;

    unsigned have_avx512;

    void *v8_buffers;           // for deallocation
    uint64_t *v8_N;             // v8_Nl elements, aligned on v8
    uint64_t *v8_a1;            // v8_Nl*2 + 8 elements, aligned on v8
    uint64_t *v8_a2;            // v8_Nl*2 + 8 elements, aligned on v8
    uint64_t v8_Nprime;         // value x: x * N mod W = -1
} Sqr_State;

// Return:
// 1: if we have needed AVX512 functionality
// 0: of no AVX512 hardware is present
static unsigned have_avx512() {
    unsigned max_function;      // largest value of eax
    unsigned eax, ebx, ecx, edx;

    __cpuid(0, eax, ebx, ecx, edx);
    max_function = eax;      // largest value of eax

    if( max_function >= 7 ) {
        SQR_TRACE("Support func 7");

        // ecx == 0 is required on input by CPUID
        __cpuid_count(7, 0, eax, ebx, ecx, edx);

        SQR_TRACE("Got ebx=%08x", ebx);

        if((ebx & ((1 << 16) | (1 << 21))) == ((1 << 16) | (1 << 21)))
            return 1;
    }
    SQR_TRACE("Return don't have AVX512");
    return 0;
}

// return a': a' * a == -1 mod W
// via Newton method. See https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
static uint64_t modinvW(uint64_t a) {
    uint64_t x = a;      //  3 bits
    x *= (2 - a * x);    //  6
    x *= (2 - a * x);    // 12
    x *= (2 - a * x);    // 24
    x *= (2 - a * x);    // 48
    x *= (2 - a * x);    // 96, capped at 64
    return (v8_W - (x & (v8_W - 1)));
}

// Transfer LE integer a of l bytes to a vector of W-sized elements out
// We assume that l is even to 8 and out has sufficient size ceil(l*8/wbits) of uint64_t elements
//
// (We don't read beyond a[l-1], as expected).
static void v8_scatter(const uint8_t *a, unsigned l, uint64_t *out) {
    unsigned i;         // input/output bit position
    uint64_t *p = (uint64_t *)a;
    // unsigned out_l = ((l * 8 + v8_wbits - 1) / v8_wbits);       // in uint64_t units
    unsigned out_l = v8_Nl;       // in uint64_t units; cover entire storage

    SQR_TRACE_BUF(a, l / 8, "Input ");

    memset(out, 0, out_l * sizeof(uint64_t));
    for(i = 0; i < l * 8; ) {
        uint64_t t = (p[i / 64] >> (i % 64)) & (v8_W - 1);
        //SQR_TRACE("p[%02d] %016jx", i/v8_wbits,  ((uint64_t)p[i / 64] >> (i % 64)));
        //SQR_TRACE("t[%02d] %016jx", i/v8_wbits, t);
        unsigned available = 64 - i % 64;
        unsigned taken = (v8_wbits - i % v8_wbits);
        if ( taken > available )
            taken = available;
        out[i / v8_wbits] |= (t << i % v8_wbits) & (v8_W - 1);
        //SQR_TRACE("out[%02d] %016jx, taken=%d (i=%d)", i / v8_wbits, out[i / v8_wbits], taken, i);

        assert(i / 64 < l / 8 && i / v8_wbits < out_l);

        i += taken;
    }
    //SQR_TRACE_BUF(out, out_l, "Return (wbits=%d, l=%d)", v8_wbits, ((l * 8 + v8_wbits - 1) / v8_wbits));
#ifdef USE_SQR_TRACE
    assert(out[((l * 8 + v8_wbits - 1) / v8_wbits)] == 0);
#endif
}

/* Pack l elements of 'a' into little-endian representation in 'out'.
 * out will be padded to (l*wbits/64)*8 bytes, which value is expected to match the input size
 * to v8_scatter.
 *
 * l is in uint64_t units
 *
 * Currently we only expect to see/tested the output size of 256 bytes (2048 bits).
 */
static void v8_gather(const uint64_t *a, unsigned l, uint8_t *out) {
    unsigned i;         // input/output bit position
    unsigned out_l = ((l * v8_wbits) / 64);       // in uint64_t units; we assume that rounding down is appropriate
    uint64_t *outp = (uint64_t *)out;

    SQR_TRACE_BUF(a, l, "Input");

    memset(outp, 0, out_l * sizeof(uint64_t));

    for(i = 0; i < out_l * 64; ) {
        uint64_t t = (a[i / v8_wbits] >> (i % v8_wbits));
        unsigned available = v8_wbits - i % v8_wbits;   // bits in t
        unsigned taken = (64 - i % 64);
        if ( taken > available )
            taken = available;
        outp[i / 64] |= (t << i % 64);
        //SQR_TRACE("out[%02d] %016jx, taken=%d (i=%d)", i / 64,  outp[i / 64], taken, i);

        assert(i / v8_wbits < l && i / 64 < out_l);

        i += taken;
    }
    SQR_TRACE_BUF(outp, out_l, "Return ");
}

/*
 * adds up to 63 bytes to p and returns this pointer, so that the returned pointer is 64-bytes aligned.
 */
static void *v8_align64_forward(void *p) {
    uintptr_t t = (uintptr_t)p;

    if((t & (64 - 1))) {
        t += 64 - (t & (64 - 1));
    }

    return (void *)t;
}

SQR_ERROR sqr_allocate_state(const uint8_t *N, unsigned Nl, Sqr_State **state) {
    BIGNUM *n = BN_new();
    BN_MONT_CTX *mont = NULL;
    BN_CTX *ctx = BN_CTX_new();
    SQR_ERROR err = SQR_ERR_OK;
    uint64_t *v8_buffers = NULL;

    *state = malloc(sizeof(Sqr_State));

    if( n == NULL || ctx == NULL || *state == NULL) {
        SQR_TRACE("Failed to allocate state or contexts");
        err = SQR_ERR_OUT_OF_MEM;
        goto error;
    }

    (*state)->have_avx512 = have_avx512();

    SQR_TRACE("Have AVX512: %d", (*state)->have_avx512);

    if( BN_bin2bn(N, Nl, n) == NULL ) {
        SQR_TRACE("Failed to import the N");
        err = SQR_ERR_BAD_PARAMS;
        goto error;
    }

    if( BN_num_bits(n) == 2048 && (*state)->have_avx512 ) {
        SQR_TRACE("Enter AVX512 code path for 2048 bit N");

        v8_buffers = malloc(v8_Nl * sizeof(uint64_t) + 2 * (2 * v8_Nl + 8) * sizeof(uint64_t) + 8 * sizeof(uint64_t));
        if( v8_buffers == NULL ) {
            SQR_TRACE("Failed to allocate v8 buffers");
            err = SQR_ERR_OUT_OF_MEM;
            goto error;
        }

        // Prepare N'
    }
    else
    {
        // default code path

        mont = BN_MONT_CTX_new();
        if( mont == NULL ) {
            err = SQR_ERR_OUT_OF_MEM;
            goto error;
        }

        if (!BN_MONT_CTX_set(mont, n, ctx)) {
            SQR_TRACE("Failed to set Montgomery context for N");
            err = SQR_ERR_INTERNAL;
            goto error;
        }
    }

    if((*state)->have_avx512 ) {
        (*state)->v8_buffers = v8_buffers;

        (*state)->v8_N = v8_align64_forward(v8_buffers);        // size Nl elements

        // working buffers of size 2*Nl + 7 elements (we use +8 instead of +7 to maintain alignment)
        (*state)->v8_a1 = (*state)->v8_N + v8_Nl;               // size 8 + 2*Nl elements
        (*state)->v8_a2 = (*state)->v8_a1 + 2 * v8_Nl + 8;      // size 8 + 2*Nl elements

        if( BN_bn2lebinpad(n, (void*)(*state)->v8_a1, 2048 / 8) < 0 ) {   // temporary the LE version of N
            err = SQR_ERR_INTERNAL;
            goto error;
        }

        // Prepare v8_N, needed for AVX512 squaring
        v8_scatter((void*)(*state)->v8_a1, 2048 / 8, (*state)->v8_N);
        //v8_gather((*state)->v8_N, (2048 + v8_wbits - 1) / v8_wbits, (void*)(*state)->v8_a2); // test

        // Calculate and store N': N' * N mod W = -1
        (*state)->v8_Nprime = modinvW((*state)->v8_N[0]);

        // lower Nl and upper 7 elements of a1 and a2 must be zero, and they will remain zero
        memset((*state)->v8_a1, 0, 2 * (2 * v8_Nl + 8) * sizeof(uint64_t));

        // avx512 success

        v8_buffers = NULL;
    }

    // complete success

    (*state)->N = n;
    (*state)->Nl = Nl;
    (*state)->montN = mont;
    (*state)->bnctx = ctx;

    n = NULL;
    mont = NULL;
    ctx = NULL;

    SQR_TRACE("Initialization OK");

    err = SQR_ERR_OK;

error:
    BN_free(n);
    BN_MONT_CTX_free(mont);
    BN_CTX_free(ctx);
    if( v8_buffers )
        free(v8_buffers);

    if( err != SQR_ERR_OK ) {
        free(*state);
        *state = NULL;
    }
    return err;
}

SQR_ERROR sqr_free_state(Sqr_State *state) {
    if( state == NULL ) {
        SQR_TRACE("Bad parameters");
        return SQR_ERR_BAD_PARAMS;
    }

    BN_free(state->N);
    BN_MONT_CTX_free(state->montN);
    BN_CTX_free(state->bnctx);
    if( state->v8_buffers )
        free(state->v8_buffers);

    free(state);

    return SQR_ERR_OK;
}


// out = x^2^T mod N
SQR_ERROR sqr_calculate(Sqr_State *state, uint64_t T, const char *xhex, char *outhex, unsigned outl) {
    BIGNUM *x = NULL;
    BIGNUM *xmont;
    unsigned l;
    unsigned i, j;
    SQR_ERROR err;

    SQR_TRACE("Begin squaring for %s, T=%08jx", xhex, T);

    if( state == NULL || T == 0 || xhex == 0 || outhex == NULL || outl == 0 || outl < state->Nl * 2 + 1) {
        SQR_TRACE("Bad parameters");
        return SQR_ERR_BAD_PARAMS;
    }

    BN_CTX_start(state->bnctx);

    x = BN_CTX_get(state->bnctx);
    if( x == NULL ) {
        SQR_TRACE("Failed to import %s", xhex);
        err = SQR_ERR_OUT_OF_MEM;
        goto error;
    }

    if( BN_hex2bn(&x, xhex) == 0 ) {
        SQR_TRACE("Failed to import %s", xhex);
        err = SQR_ERR_BAD_PARAMS;
        goto error;
    }

    l = BN_num_bytes(x);
    if( l > state->Nl ) {
        SQR_TRACE("Internal consistency error: %d > %d", l, state->Nl);
        err = SQR_ERR_INTERNAL;
        goto error;
    }

    xmont = BN_CTX_get(state->bnctx);
    if( xmont == NULL ) {
        SQR_TRACE("Failed to allocate a BN");
        err = SQR_ERR_OUT_OF_MEM;
        goto error;
    }

    if(!state->have_avx512 ) {

        // Default code path
        if( !BN_to_montgomery(xmont, x, state->montN, state->bnctx)) {
            SQR_TRACE("Failed to convert to Montgomery");
            err = SQR_ERR_BAD_PARAMS;
            goto error;
        }

        while( T-- )  {
            if (!BN_mod_mul_montgomery(xmont, xmont, xmont, state->montN, state->bnctx)) {
                err = SQR_ERR_OUT_OF_MEM;
                goto error;
            }
        }

        if (!BN_from_montgomery(x, xmont, state->montN, state->bnctx)) {
            err = SQR_ERR_OUT_OF_MEM;
            goto error;
        }

        // write in binary, padded to Nl, at the end of the output buffer outhex
        j = outl - state->Nl; // start of binary x
        if( BN_bn2binpad(x, (unsigned char *)outhex + outl - state->Nl, state->Nl) <= 0 ) {
            SQR_TRACE("Internal consistency error");
            err = SQR_ERR_INTERNAL;
            goto error;
        }

        // convert to hex in-place
        for( i = 0; i < state->Nl; i++ ) {
            snprintf(outhex + i * 2, 3, "%02x", (uint8_t)outhex[j++]);
        }

        assert(j == outl);
    }
    else
    {
        // AVX512 code path for 2048-bit N

        SQR_TRACE("In AVX512 code path");

        // To our Montgomery form.
        // Calculate xmont = x * W^(l/8) mod N with an undocumented function, which
        // is perfect for us.
        // BN_mod_lshift_quick requires 0 < x < N, and does combined left shift and reduction
        if( !BN_mod_lshift_quick(xmont, x, v8_wbits * v8_Nl_tight, state->N)) {
            SQR_TRACE("Failed to convert to Montgomery");
            err = SQR_ERR_INTERNAL;
            goto error;
        }
        SQR_TRACE_BN(xmont, "xmont for %s", xhex);

        // convert xmont to expanded array
        if( BN_bn2lebinpad(xmont, (uint8_t*)(state->v8_a2 + v8_Nl), state->Nl) <= 0 ) {       // xmont --> v8_a2, packed
            SQR_TRACE("Failed to serialize a BN");
            err = SQR_ERR_INTERNAL;
            goto error;
        }
        SQR_TRACE_BUF(state->v8_a2 + v8_Nl,  state->Nl / 8, "xmont packed input (a2)");
        v8_scatter((const uint8_t*)(state->v8_a2 + v8_Nl), 2048 / 8, state->v8_a1 + v8_Nl);  // a2 --> a1 (upper Nl elements)
        SQR_TRACE_BUF(state->v8_a1, 2 * v8_Nl, "xmont expanded into second Nl elements (a1)");

        SQR_TRACE("Begin squaring loop, 2^%08jx", T);

        v8_square_2048(T, state->v8_a1 + v8_Nl, state->v8_N, state->v8_Nprime, state->v8_a2);   // a1^2 --> a2 (second Nl elements)

        SQR_TRACE_BUF(state->v8_a2, 2 * v8_Nl + 8, "xmont after loop");
        v8_from_mont(state->v8_a2 + v8_Nl, state->v8_N, state->v8_Nprime, v8_Nl_tight, state->v8_a1 + v8_Nl);         // mont_a1 --> a2
        SQR_TRACE_BUF(state->v8_a1 + v8_Nl, v8_Nl, "x after loop");

        v8_gather(state->v8_a1 + v8_Nl, v8_Nl, (uint8_t*)(state->v8_a2 + v8_Nl));  // a2 --> a1, packed LE
        SQR_TRACE_BUF(state->v8_a2, 2 * v8_Nl, "x after loop, packed");

        const uint8_t * p = (const uint8_t*)(state->v8_a2 + v8_Nl);

        // convert from LE a1 to hex BE in outhex
        for( i = 0; i < 2048 / 8; i++ ) {
            snprintf(outhex + i * 2, 3, "%02x", p[2048 / 8 - 1 - i]);
        }
    }

    // Success
    err = SQR_ERR_OK;

error:
    BN_CTX_end(state->bnctx);

    return err;
}


#include <stdio.h>
#include <openssl/bn.h>
#include <assert.h>

#include "sqrlib.h"

#if 1
#define SQR_TRACE(format, ...) \
    printf("%s(%d): " format "\n", __FILE__, __LINE__, ## __VA_ARGS__)
#else
#define SQR_TRACE(format, ...)
#endif


typedef struct Sqr_State {
    BIGNUM *N;
    unsigned Nl;
    BN_MONT_CTX *montN;
    BN_CTX *bnctx;
} Sqr_State;

SQR_ERROR sqr_allocate_state(const uint8_t *N, unsigned Nl, Sqr_State **state) {
    BIGNUM *n = BN_new();
    BN_MONT_CTX *mont = BN_MONT_CTX_new();
    BN_CTX *ctx = BN_CTX_new();
    SQR_ERROR err = SQR_ERR_OK;

    *state = malloc(sizeof(Sqr_State));

    if( n == NULL || mont == NULL || ctx == NULL || *state == NULL) {
        SQR_TRACE("Failed to allocate state or contexts");
        err = SQR_ERR_OUT_OF_MEM;
        goto error;
    }
    if( BN_bin2bn(N, Nl, n) == NULL ) {
        SQR_TRACE("Failed to import the N");
        err = SQR_ERR_BAD_PARAMS;
        goto error;
    }

    if (!BN_MONT_CTX_set(mont, n, ctx)) {
        SQR_TRACE("Failed to set Montgomery context for N");
        goto error;
    }

    // success

    (*state)->N = n;
    (*state)->Nl = Nl;
    (*state)->montN = mont;
    (*state)->bnctx = ctx;

    n = NULL;
    mont = NULL;
    ctx = NULL;

    err = SQR_ERR_OK;

error:
    BN_free(n);
    BN_MONT_CTX_free(mont);
    BN_CTX_free(ctx);

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

    if( state == NULL || T == 0 || xhex == 0 || outhex == NULL || outl == 0 || outl < state->Nl * 2 + 1) {
        SQR_TRACE("Bad parameters");
        return SQR_ERR_BAD_PARAMS;
    }

    BN_CTX_start(state->bnctx);

    x = BN_CTX_get(state->bnctx);

    if( BN_hex2bn(&x, xhex) == 0 ) {
        SQR_TRACE("Failed to import %s", xhex);
        err = SQR_ERR_BAD_PARAMS;
        goto error;
    }

    xmont = BN_CTX_get(state->bnctx);

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

    l = BN_num_bytes(x);
    if( l > state->Nl ) {
        SQR_TRACE("Internal consistency error: %d > %d", l, state->Nl);
        err = SQR_ERR_INTERNAL;
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

    // Success
    err = SQR_ERR_OK;

error:
    BN_CTX_end(state->bnctx);

    return err;
}


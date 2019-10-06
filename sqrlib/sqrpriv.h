#ifndef SQR_PRIVATE_H
#define SQR_PRIVATE_H 1

#define v8_Nl 48                // length in elements; element = uint64_t word
#define v8_wbits 50             // actual bits in the element, the remaining are for carries; 50*41 = 2050 (41 actual elements: 7 high elements in the last v8 are zeros)
#define v8_W (1ULL << v8_wbits)
#define v8_Nl_tight 41       // minimum such that > 2048

void v8_square_2048(uint64_t T, const uint64_t a[v8_Nl], const uint64_t N[v8_Nl], uint64_t Nprime, uint64_t out[7+2*v8_Nl]);

void v8_from_mont(const uint64_t *a, const uint64_t *N, uint64_t Nprime, unsigned l, uint64_t *a_out);

//#define USE_SQR_TRACE 1

#ifdef USE_SQR_TRACE
#define SQR_TRACE(format, ...) \
    printf(". %s(%d):%s: " format "\n", __FILE__, __LINE__, __func__, ## __VA_ARGS__)

static inline void print(const void *v, unsigned l) {
    unsigned i;
    const uint64_t *p = (const uint64_t *)v;
    for( i = 0; i < l; i++ ) {
        if( i % 8 == 0 ) {
            printf(".    ");
        }
        printf("%016jx", p[i]);
        if( i != l - 1 )
            printf(", ");
        if( i % 8 == 7 && i != l - 1 ) {
            printf("\n");
        }
    }
    printf("\n");
}

#define SQR_TRACE_BUF(vector, l, format, ...) \
    printf(". %s(%d):%s: " format " [%d]\n", __FILE__, __LINE__, __func__, ## __VA_ARGS__, l); \
    print(vector, l);

#define SQR_TRACE_BN(bn, format, ...) \
    { \
        uint8_t out[2048 / 8]; \
        if( BN_bn2lebinpad(bn, out, sizeof(out)) < 0 )  /* LE version of bn */ \
        printf(". %s(%d):%s: " format "FAILED\n", __FILE__, __LINE__, __func__, ## __VA_ARGS__); \
        else { \
            SQR_TRACE_BUF((void*)out, (unsigned)sizeof(out) / 8, format,  ## __VA_ARGS__); \
        } \
    }

#else
#define SQR_TRACE(format, ...)
#define SQR_TRACE_BUF(buf, format, ...)
#define SQR_TRACE_BN(bn, format, ...)
#endif

#endif // SQR_PRIVATE_H

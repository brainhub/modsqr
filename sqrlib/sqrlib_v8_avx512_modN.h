#define SHR1_N_6() \
    N8_5 = _mm512_alignr_epi64(N8_5, N8_4, 8 - 1); \
    N8_4 = _mm512_alignr_epi64(N8_4, N8_3, 8 - 1); \
    N8_3 = _mm512_alignr_epi64(N8_3, N8_2, 8 - 1); \
    N8_2 = _mm512_alignr_epi64(N8_2, N8_1, 8 - 1); \
    N8_1 = _mm512_alignr_epi64(N8_1, N8_0, 8 - 1); \
    N8_0 = _mm512_alignr_epi64(N8_0, zero, 8 - 1);

// a*b --> {l, h}: h*2^52 + l
#define MULX(a,b,mask,outl, outh) \
    asm ("mulx %3, %0, %1; " \
         "shld $(64-52), %0, %1; " \
         "and %4, %0;" \
         "and %4, %1;" \
         : "=r" (outl), "=r" (outh) /* output #0 for low, #1 for high*/ \
         : "d" (a), "r" (b), "r" (mask) /* input #2 is in %rdx, input #3 is any register */ \
         : "cc" /* clobbered registers: flags */ \
         );

#define RED_N(i, _O_0,_O_1,_O_2,_O_3,_O_4,_O_5,_O_6,  _Oh_0,_Oh_1,_Oh_2,_Oh_3,_Oh_4,_Oh_5) \
    { \
        __m512i M8; \
        const uint64_t n0 = N[0]; /*faster than _mm_extract_epi64(_mm512_extracti32x4_epi32(N8_0, 0), 0);*/ \
        const uint64_t n40 = _mm_extract_epi64(_mm512_extracti32x4_epi32(N8_5, 0), 0); /* faster than N[40] */ \
        \
        uint64_t t0 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 0), 0); /* O_0[0]*/ \
        uint64_t th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 0), 0) << (52 - wbits); \
        uint64_t t1 = t0 + th1; \
        uint64_t m = (t1 * Nprime) & (W - 1); \
        const uint64_t mask = (1ULL << 52) - 1; \
        uint64_t l, h; \
        \
        M8 = _mm512_set1_epi64(m); \
\
        /* j=0  */ \
        /*_O_0 = _mm512_mask_add_epi64(_O_0, 1, _O_0, _mm512_set1_epi64(th1));*/ \
\
        MULX(n0, m, mask, l, h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 0), 1); /* Oh_0[1] */ \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
        \
        SHR1_N_6();       /* N<<=64 */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 0), 1); /* O_0[1] */ \
\
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0);   \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
\
        /* j=1 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 1), 0); \
\
        SHR1_N_6();       /* N<<=64 (2*64 at this point) */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 1), 0); \
\
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
\
        /* j=2 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 1), 1); \
        \
        SHR1_N_6(); \
        \
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 1), 1); \
\
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
\
        /* j=3 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 2), 0); \
        \
        SHR1_N_6();       /* N<<=64 */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 2), 0); \
        \
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
\
        /* j=4 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 2), 1); \
        \
        SHR1_N_6();       /* N<<=64 */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 2), 1); \
        \
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
       \
        /* j=5 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
 \
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 3), 0); \
        \
        SHR1_N_6();       /* N<<=64 */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 3), 0); \
        \
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, N8_0); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
        \
        /* j=6 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        MULX(n0,m,mask,l,h); \
        t0 = (t1 + l) >> wbits; \
        th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 3), 1); \
        \
        SHR1_N_6();       /* N<<=64 */ \
\
        t1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 3), 1); \
        \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_1); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_2); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_3); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_4); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_5); \
\
        /* j=7 */ \
        t1 = (t1 + t0 + ((th1 + h) << (52 - wbits))); m = (t1 * Nprime) & (W - 1); M8 = _mm512_set1_epi64(m); \
\
        asm ("mulx %3, %1, %0; " \
             "shrd $52, %0, %1; " /* wbits=50 */ \
             "shl $2, %1;" \
             : "=r" (l), "=r" (h) /* output #0 for low (dummy), #1 for high*/ \
             : "d" (m), "r" (n40) /* input #2 is in %rdx, input #3 is any register */ \
             : "cc" /* clobbered registers: flags */ \
             ); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
\
        /* restore N */ \
        N8_0 = _mm512_alignr_epi64(N8_1, N8_0, 7);  \
        N8_1 = _mm512_alignr_epi64(N8_2, N8_1, 7); \
        N8_2 = _mm512_alignr_epi64(N8_3, N8_2, 7); \
        N8_3 = _mm512_alignr_epi64(N8_4, N8_3, 7); \
        N8_4 = _mm512_alignr_epi64(N8_5, N8_4, 7); \
        N8_5 = _mm512_maskz_set1_epi64(1, n40); \
        \
        t0 = (t1 + ((n0 * m) & mask)) >> wbits; \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, N8_0);   /* full shift */ \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, N8_1); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, N8_2); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, N8_3); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, N8_4); \
        \
        _O_6 = _mm512_mask_add_epi64(_O_6, 1, _O_6, _mm512_set1_epi64(h)); \
       \
        _O_1 = _mm512_mask_add_epi64(_O_1, 1, _O_1, _mm512_set1_epi64(t0)); \
    }

/* make O_0[0] == 0 by adding m*N once */
#define RED_N_1(i, _O_0,_O_1,_O_2,_O_3,_O_4,_O_5,  _Oh_0,_Oh_1,_Oh_2,_Oh_3,_Oh_4,_Oh_5) \
    { \
        __m512i M8; \
        uint64_t t0 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_O_0, 0), 0); /* O_0[0]*/ \
        uint64_t th1 = _mm_extract_epi64(_mm512_extracti32x4_epi32(_Oh_0, 0), 0) << (52 - wbits); \
        uint64_t t1 = t0 + th1; \
        uint64_t m = (t1 * Nprime) & (W - 1); \
        const uint64_t mask = (1ULL << 52) - 1; \
        \
        M8 = _mm512_set1_epi64(m); \
\
        /* j=0  */ \
\
        t0 = (t1 + ((N[0] * m) & mask)) >> wbits; \
\
        _O_0 = _mm512_madd52lo_epu64(_O_0, M8, N8_0); \
        _O_1 = _mm512_madd52lo_epu64(_O_1, M8, N8_1); \
        _O_2 = _mm512_madd52lo_epu64(_O_2, M8, N8_2); \
        _O_3 = _mm512_madd52lo_epu64(_O_3, M8, N8_3); \
        _O_4 = _mm512_madd52lo_epu64(_O_4, M8, N8_4); \
        _O_5 = _mm512_madd52lo_epu64(_O_5, M8, N8_5); \
        \
        _Oh_0 = _mm512_madd52hi_epu64(_Oh_0, M8, _mm512_alignr_epi64(N8_0, zero, 8 - 1)); \
        _Oh_1 = _mm512_madd52hi_epu64(_Oh_1, M8, _mm512_alignr_epi64(N8_1, N8_0, 8 - 1)); \
        _Oh_2 = _mm512_madd52hi_epu64(_Oh_2, M8, _mm512_alignr_epi64(N8_2, N8_1, 8 - 1)); \
        _Oh_3 = _mm512_madd52hi_epu64(_Oh_3, M8, _mm512_alignr_epi64(N8_3, N8_2, 8 - 1)); \
        _Oh_4 = _mm512_madd52hi_epu64(_Oh_4, M8, _mm512_alignr_epi64(N8_4, N8_3, 8 - 1)); \
        _Oh_5 = _mm512_madd52hi_epu64(_Oh_5, M8, _mm512_alignr_epi64(N8_5, N8_4, 8 - 1)); \
        \
        _O_0 = _mm512_mask_add_epi64(_O_0, 1 << 1, _O_0, _mm512_set1_epi64(t0)); \
    }




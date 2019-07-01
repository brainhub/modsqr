#include <cpuid.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef uint8_t u8;

#define UtoC(u) (u8)u,(u8)(u >> 8),(u8)(u >> 16),u >> 24

void print_cpuid_info() {
    unsigned max_function;      // largest value of eax
    unsigned eax, ebx, ecx, edx;
    unsigned have_sse2 = 0;
    unsigned have_avx = 0;
    unsigned have_mulx_adx = 0;

    __cpuid(0, eax, ebx, ecx, edx);
    max_function = eax;      // largest value of eax

    printf("CPU information:\n");
    printf("  %-16s: %c%c%c%c" "%c%c%c%c" "%c%c%c%c\n",
           "Brand",
           UtoC(ebx), UtoC(edx), UtoC(ecx)
           );

    __cpuid(0x80000000, eax, ebx, ecx, edx);
    if( eax >= 0x80000004) {
        char s[49];
        char *p = s;

        __cpuid(0x80000002, eax, ebx, ecx, edx);
        snprintf(s, 17,
                 "%c%c%c%c" "%c%c%c%c" "%c%c%c%c" "%c%c%c%c",
                 UtoC(eax), UtoC(ebx), UtoC(ecx), UtoC(edx));

        __cpuid(0x80000003, eax, ebx, ecx, edx);
        snprintf(s + 16, 17,
                 "%c%c%c%c" "%c%c%c%c" "%c%c%c%c" "%c%c%c%c",
                 UtoC(eax), UtoC(ebx), UtoC(ecx), UtoC(edx));

        __cpuid(0x80000004, eax, ebx, ecx, edx);
        snprintf(s + 32, 17,
                 "%c%c%c%c" "%c%c%c%c" "%c%c%c%c" "%c%c%c%c",
                 UtoC(eax), UtoC(ebx), UtoC(ecx), UtoC(edx));

        while(*p == ' ')
            p++;

        printf("  %-16s: %s\n", "Brand", p);
    }

    //printf( "Basic max function: %d\n", max_function );

    if( max_function >= 1 ) {
        __cpuid(1, eax, ebx, ecx, edx);
        //  OPENSSL_ia32cap_P[0]=edx  OPENSSL_ia32cap_P[1]=ecx for func #1
        printf("  %-16s: ebx=%08x ecx=%08x edx=%08x\n", "Features", ebx, ecx, edx);

        // Older top-level flags are in EDX, func #1
        // ECX adds a few more, as for AVX and AES.
        // They are in OPENSSL_ia32cap_P[0]=edx, OPENSSL_ia32cap_P[1]=ecx and controlled by e.g.
        // OPENSSL_ia32cap=~0x4000000 apps/openssl speed rsa3072
        // to disable SSE2. In my tests disabling SSE2 doesn't matter with 1.1.0 code.
        //
        // Same for AVX. AVX can be disabled as
        // OPENSSL_ia32cap=~0x1000000000000000  apps/openssl speed rsa3072
        // This also doesn't matter for this code.
        //
        // For reference, the following disables AESNI, which location is close to AVX:
        // OPENSSL_ia32cap="~0x100000000000000" openssl speed -elapsed -evp aes-128-cbc
        //
        // Stop at SSE2 and AVX for now. You should always see them, and they indicate
        // baseline support on all systems
        have_sse2 = (edx & 0x4000000) == 0x4000000;     // (1<<26)
        have_avx = (ecx & 0x10000000) == 0x10000000;    // (1<<28)
    }

    if( max_function >= 7 ) {
        // ecx == 0 is also required on input
        __cpuid(7, eax, ebx, ecx, edx);
        // OPENSSL_ia32cap_P[2]=281=ebx  OPENSSL_ia32cap_P[3]=0=ecx for func #7

        printf("  %-16s: ebx=%08x ecx=%08x edx=%08x\n", "Ext. features", ebx, ecx, edx);

        // Disable MULX, ADX: OPENSSL_ia32cap=":~0x80100"  apps/openssl speed rsa3072
        //
        // Force-enable:
        // OPENSSL_ia32cap=":0x80100"  apps/openssl speed rsa3072
        // or this for Ivy Bridge
        // OPENSSL_ia32cap=":0x80381"  apps/openssl speed rsa3072
        // This will crash this code, as expected.

        // ebx 0x80100:
        // 0x100 : MULX
        // 0x80000 : ADX
        have_mulx_adx = (ebx & 0x80100) == 0x80100;
    }
    printf("  %-16s: %s\n", "SSE2", have_sse2 ? "yes" : "no");
    printf("  %-16s: %s\n", "AVX", have_avx ? "yes" : "no");

    printf("  %-16s: %s\n", "MULX, ADX", have_mulx_adx ? "yes" : "no (IMPORTANT!)");

    putchar('\n');
}

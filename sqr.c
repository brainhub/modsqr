#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "sqrlib.h"

// N=25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357
// from https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048
static const uint8_t RSA_2048_N[] = {
    0xc7,0x97,0x0c,0xee, 0xdc,0xc3,0xb0,0x75, 0x44,0x90,0x20,0x1a, 0x7a,0xa6,0x13,0xcd, 0x73,0x91,0x10,0x81, 0xc7,0x90,0xf5,0xf1, 0xa8,0x72,0x6f,0x46, 0x35,0x50,0xbb,0x5b,
    0x7f,0xf0,0xdb,0x8e, 0x1e,0xa1,0x18,0x9e, 0xc7,0x2f,0x93,0xd1, 0x65,0x00,0x11,0xbd, 0x72,0x1a,0xee,0xac, 0xc2,0xac,0xde,0x32, 0xa0,0x41,0x07,0xf0, 0x64,0x8c,0x28,0x13,
    0xa3,0x1f,0x5b,0x0b, 0x77,0x65,0xff,0x8b, 0x44,0xb4,0xb6,0xff, 0xc9,0x33,0x84,0xb6, 0x46,0xeb,0x09,0xc7, 0xcf,0x5e,0x85,0x92, 0xd4,0x0e,0xa3,0x3c, 0x80,0x03,0x9f,0x35,
    0xb4,0xf1,0x4a,0x04, 0xb5,0x1f,0x7b,0xfd, 0x78,0x1b,0xe4, 0xd1,0x67,0x31,0x64,0xba, 0x8e,0xb9,0x91, 0xc2,0xc4,0xd7,0x30,0xbb, 0xbe,0x35,0xf5,0x92, 0xbd,0xef,0x52,0x4a,
    0xf7,0xe8,0xda,0xef, 0xd2,0x6c,0x66,0xfc, 0x02,0xc4,0x79,0xaf, 0x89,0xd6,0x4d,0x37, 0x3f,0x44,0x27,0x09, 0x43,0x9d,0xe6,0x6c, 0xeb,0x95,0x5f,0x3e, 0xa3,0x7d,0x51,0x59,
    0xf6,0x13,0x58,0x09, 0xf8,0x53,0x34,0xb5, 0xcb,0x18,0x13,0xad, 0xdc,0x80,0xcd,0x05, 0x60,0x9f,0x10,0xac, 0x6a,0x95,0xad,0x65, 0x87,0x2c,0x90,0x95, 0x25,0xbd,0xad,0x32,
    0xbc,0x72,0x95,0x92, 0x64,0x29,0x20,0xf2, 0x4c,0x61,0xdc,0x5b, 0x3c,0x3b,0x79,0x23, 0xe5,0x6b,0x16,0xa4, 0xd9,0xd3,0x73,0xd8, 0x72,0x1f,0x24,0xa3, 0xfc,0x0f,0x1b,0x31,
    0x31,0xf5,0x56,0x15, 0x17,0x28,0x66,0xbc, 0xcc,0x30,0xf9,0x50, 0x54,0xc8,0x24,0xe7, 0x33,0xa5,0xeb,0x68, 0x17,0xf7,0xbc,0x16, 0x39,0x9d,0x48,0xc6, 0x36,0x1c,0xc7,0xe5
};

void static print_help() {
        printf("Square x 2^t times. To call:\n"
            "sqr <t: decimal integer> <x: hexidecimal integer>\n"
            "\n"
            "Example:\n"
            "\n"
            "$ sqr 40 8fc21da610761c189aa419477fef0211000b04583f10abda004ee47d2e98626c\n"
            ""
        );
}

int main(int argc, const char * argv[])
{
    Sqr_State *state;
    char out[2 * 2048 / 8 + 1];
    unsigned t;
    unsigned xlen;
    SQR_ERROR err;

    if( argc != 3 ) {
        print_help();
        return 1;
    }

    t = atoi(argv[1]);
    if( t <= 0 || t>63 ) {
        printf("t is too large. 63 is maximum\n");
        return 1;
    }

    xlen = strlen(argv[2]);
    
    if( xlen < 256/8 || xlen > 2048/8*2 ) {
        if( xlen < 256/8 )
            printf("The x is too short (smaller than 128-bit integer)\n");
        if( xlen > 2048/8*2 )
            printf("The x is too large (larger than 2048-bit integer)\n");
        return 1;
    }

    err = sqr_allocate_state(RSA_2048_N, sizeof(RSA_2048_N), &state);
    if( err != SQR_ERR_OK )
        return 1;

    err = sqr_calculate(state, (uint64_t)1ULL<<t, argv[2], out, sizeof(out));
    sqr_free_state(state);
    if( err != SQR_ERR_OK ) {
        printf("Failed to calculate x^2^(2^t)\n");
        return 1;
    }

    puts(out);

    return 0;
}
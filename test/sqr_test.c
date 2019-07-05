#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include <x86intrin.h>

// for version info only
#include <openssl/opensslv.h>
#include <openssl/crypto.h>

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

const char * const mhex = "5bc382da0c21ce5dbac814bd382ac3bee474e467500800e32cbd7f783adf0940";  // sha256 sum of "vdf"

// (m^2) mod N
// Observe mandatory padding to N bytes (fixed size text output)
const char * const m_2 =
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000020e494579fb0b2420552a4b3140d5c1157f1d55c9506d9dfe161228dcf7d8e4b0c49a3ba23e0c58ff81163f977b09e5286244487c6c682e2324ffe821dd59000";

// (m^2) mod N
const char * const m_2_2 =
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000439f12e79f60df042c09accea1dea5dfd5031a690b773d31e49b8bce0109673a2d1216b9ceae9c909ca5e5903f51a410fc2d43385cae66b489950da94902063360f79e51b200729da14f9bc227785d450488f0b38edeeedb3b046bbcddaafab052b08ebc8e917a28c672054b49bd84a80c73d8e384e8b64fa2054c8f1000000";

// (m^(2^40)) mod N
const char * const m_2_40 =
    "9d0073ae84cd3517d3db2f2b0000a82737248b576d8df421549d19320e7823a2c706dd9bdf21c4bf940aa870fd10465b3ea1f6131d629a4a849af14ee2afce19d905379963f93d39b19768ad23445a9a6d23111fad6335c1db3f5601f2151894e3034cac5f7c6f335f7a8a425805a27a152196f4437bc5424139a0fbd22201531e32adb71db66f5b9125e127cae4d9790a2174e81f1a019f6bf21fd9430a308f6bdd0241604f2bed98a37e4a169e15b376e3eb3f0e48cefaaad1fa17d63abdcdc173b72047a65652fe35ec421200ad72fcb13f492976cf93052ed741a14a2c343974e4dd1fb600329417fbd70e49139c562dfea58af68fdd1ab399639404d1ec";

// m^(2^2^10) == m^2^1024 mod N
const char * const m_2_2_20 =
    "4ed421fa52dd1b828839b74b07f8ce6de0a16f752c893cf49dc55c23fd9ff3da5538df0979f42bc5b9a2dbb14c25b2da5a57ef814995ffba82441aac9633d537bddabd5f6229aff0cf92059aa3d22837e9daa06182968cc329f5926d53ac701d75add55494ced38da658faa7737aa19e49188b7a2621174e8897665dbe0d79c6f45b917112b5ff47b75aa45b4a8b1c0372185190e229e549fc58adfb1a452bf146bcbef3970737abd2e6a928baeb6c1b45f5f0a4188dc81e8d61f5f6413f3f5d3f03337eb8cf850e6388fd3a61fd33a831d552f6cecde6ef16f184dd91a0c717cfd81aa604277e799e46b09d0e55224503d7c9711596c88e7af73500afa30cf1";

static inline uint64_t cycles()
{
    return __rdtsc();
}

static uint64_t time_now() {
    struct timeval tv;
    uint64_t ret;

    gettimeofday(&tv, NULL);
    ret = tv.tv_sec;
    ret *= 1000000;
    ret += tv.tv_usec;

    return ret;
}

void print_cpuid_info();

// ceiling for log2(n);
static unsigned l2(uint64_t n) {
    uint64_t mask =  (uint64_t)1ULL << 63;
    unsigned r = 63;

    if( n <= 1 )
        return 1;

    while( !(mask & n))  {
        mask >>= 1;
        r--;
    }

    return r + !!(n & (n - 1));
}


int main(int argc, const char * argv[])
{
    Sqr_State *state;
    char out[2 * 2048 / 8 + 1];
    uint64_t tbegin, tend;
    uint64_t cbegin, cend;
    uint32_t estimate_time, actual_time;
    unsigned latest_t;
    SQR_ERROR err;

    print_cpuid_info();

    printf( "OpenSSL version: %s\n", OpenSSL_version(OPENSSL_VERSION));
    if( OpenSSL_version_num() < 0x10100000 ) {  // < 1.1.x ? 
        printf( "WARNING: OpenSSL version is too low (%08lx)\n", OpenSSL_version_num());
    }

    putchar('\n');

    err = sqr_allocate_state(RSA_2048_N, sizeof(RSA_2048_N), &state);
    if( err != SQR_ERR_OK )
        return 1;

    err = sqr_calculate(state, 1, mhex, out, sizeof(out));
    if( strcmp(out, m_2) != 0 ) {
        printf("Failed in x^2:\n%s != \n%s\n", out, m_2);
        sqr_free_state(state);
        return 1;
    }

    err = sqr_calculate(state, 2, mhex, out, sizeof(out));
    if( strcmp(out, m_2_2) != 0 ) {
        printf("Failed in x^2^2:\n%s != \n%s\n", out, m_2_2);
        sqr_free_state(state);
        return 1;
    }

    err = sqr_calculate(state, 40, mhex, out, sizeof(out));
    if( strcmp(out, m_2_40) != 0 ) {
        printf("Failed in x^2^40:\n%s != \n%s\n", out, m_2_40);
        sqr_free_state(state);
        return 1;
    }

    latest_t = 20;
    tbegin = time_now();
    cbegin = cycles();
    err = sqr_calculate(state, ((uint64_t)1ULL) << latest_t, mhex, out, sizeof(out));
    if( latest_t == 20 && strcmp(out, m_2_2_20) != 0 ) {
        printf("Failed in x^2^2^%d:\n%s != \n%s\n", latest_t, out, m_2_2_20);
        sqr_free_state(state);
        return 1;
    }
    tend = time_now();
    cend = cycles();
    printf("%g op/sec in %.02f sec for x^2^2^%d\n",  (1 << latest_t) * 1000000. / (tend - tbegin), (tend - tbegin) / 1000000., latest_t);
    printf("cycles for one square mod N: %lu\n", (cend - cbegin) / (1 << latest_t));

    putchar('\n');

    estimate_time = (unsigned)((tend - tbegin) / 1000000. * 32.);

    if( estimate_time < 5 * 60 ) {
        printf("Estimated time to complete the next test is %d seconds\n", estimate_time);
        tbegin = time_now();
        cbegin = cycles();
        err = sqr_calculate(state, ((uint64_t)1ULL) << 25, mhex, out, sizeof(out));
        assert(err == SQR_ERR_OK);
        tend = time_now();
        cend = cycles();
        printf("%g op/sec in %.02f sec for x^2^2^25\n",  (1 << 25) * 1000000. / (tend - tbegin), (tend - tbegin) / 1000000.);
        printf("cycles for one square mod N: %lu\n", (cend - cbegin) / (1 << 25));
        actual_time =  (unsigned)((tend - tbegin) / 1000000. * 32.);
        latest_t = 25;

        if( estimate_time > actual_time + actual_time / 10 || estimate_time > actual_time - actual_time / 10 ) {
            printf("Warning: estimated time to complete differs by over 10%% from the actual time. Don't use TurboBoost, etc. "
                   "Following estimates may be inaccurate.\n");
        }

        putchar('\n');
    }
    else
    {
        printf("Skipping x^2^2^25 calculation because it will take more then 5 min\n");
    }

    /*
     * It takes (tend - tbegin) / 1000000. sec to do 2^25 operations.
     * Calculate t needed to delay for one day.
     */
    {
        unsigned times_per_day_est = 3600. * 24. * 1000000. / (tend - tbegin);
        unsigned bits = l2(times_per_day_est);

        assert(1 << bits >= times_per_day_est);
        //printf("2^%d >= %g\n", bits, (float)times_per_day_est );

        printf("x^2^2^%d will take approximately 1 day and %g min (%g day) to compute on this system (t=%d)\n",
               latest_t + bits,
               (1 << bits) * (tend - tbegin) / (60. * 1000000.) - 24. * 60.,
               (1 << bits) * (tend - tbegin) / (24. * 3600. * 1000000.),
               latest_t + bits);
    }

    sqr_free_state(state);

    return 0;
}

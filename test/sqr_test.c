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

/* m^(2^2^20) mod N
    N=25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357
    R = Integers(N)
    m = 0x5bc382da0c21ce5dbac814bd382ac3bee474e467500800e32cbd7f783adf0940
    m1 = R(m)
    for i in range(1024*1024):
        m1 = m1*m1
    hex(int(m1))
 # matches the below m_2_2_20
 */
const char * const m_2_2_20 =
    "4ed421fa52dd1b828839b74b07f8ce6de0a16f752c893cf49dc55c23fd9ff3da5538df0979f42bc5b9a2dbb14c25b2da5a57ef814995ffba82441aac9633d537bddabd5f6229aff0cf92059aa3d22837e9daa06182968cc329f5926d53ac701d75add55494ced38da658faa7737aa19e49188b7a2621174e8897665dbe0d79c6f45b917112b5ff47b75aa45b4a8b1c0372185190e229e549fc58adfb1a452bf146bcbef3970737abd2e6a928baeb6c1b45f5f0a4188dc81e8d61f5f6413f3f5d3f03337eb8cf850e6388fd3a61fd33a831d552f6cecde6ef16f184dd91a0c717cfd81aa604277e799e46b09d0e55224503d7c9711596c88e7af73500afa30cf1";

typedef struct {
    unsigned t;
    const char *y;
} KAT;

KAT kat[] = {
    { 5, "259a3eacd1c07fe34f54fbe84f963ae4bc0aaeffff25a32f849c4163a0fbf7bcbdc527891ddd1a8a4d691f3f3c8472e33c9c2675d62dbe7aca60ab021dddada8a278f0439790bcad7eeca951110d88ec2c771ac65530a65a2a6c6a050397fde790db0fb3dc33960f56757afa728cd1568538fd848fd4a2cf8cb6cd23c77bb185d1fa3ab50ff08388e416fec081e4f12028dd044f00baf60343e0607f46ddc6ef2ececb7d2b42258d7aa775c8d50626e3d10e84dafb86e20102ee6f2ee4dc35969069c800f3f3f208ef3c539decf3aa7cd4d8c51e9b32185faaa6fa2e8a52682618e725ad380bc0ec7517426ba383ef8cfd9e45b392d8bd1d64f3048ba1cb86bd"},
    {10, "9f5cea72f596defe89f60ff4e5e6c018f4e145eba09db647a0d8f145dcf06230a257dcd0655fb8d55293f2fec26003f3dca90f6b5fbd38db4829ecd9523264279a8487cd066478a3b1cb33764ecea940aa74ee54362009fa0ad6ae75935c6c5beed7ac74fa28ce3e38d511e72b7c21b9f6eda51af6d086078006584bc1e3c5808ed4ec4bd8a4ebe48d6b537d51b80888a4fa380a42d8058394d6353db2b651d51da2efd9ff8451ec3a320cfeb746ac1860f191cbcf7bd046acd50183a545e2197b81936fd41f97ae729ce6ef9b1ac13ea6d269b89ce9eb1ddea895976f01d66d2e45859054e1514ee81ec17a7260fee59a03f83919aeaf5efb4b2e0be6b76541"},
    {11, "bd458261df583f9adf905527fba3149f150ac3f7c777cdbb034541330de0b85e9405306f8067f35c8797aee7ed46d13280e32c90a4bb70a52a1f49da8a32f5c8cbfb0e14a7375c102f26b9abb8ee517e092439140e44a3f97ae7028c4febb813e42081b00244e41f0ae1987f4d582b26499f8236a926fd3cbe955ddefaaebebca405fbd0138a3739b4e7b1de3f142233a94d732be3f23beacf9dc0b00052e48152820192c7ecbbd4cb2ed2e7b377b693bec28f636283e205591238d4c3e8e72af878fc2db0e0a07dbc7bc3d190279a630f353e89e09bcc18e6d00b33a3ae0ebc7ef9a070d53e816e0084833eb2d29b3eb780c353a698795a687830fc11c15e75"},
    {12, "8870f1732ec0e6f233b5f6883b804bd5fc588b384b5baad872741b2917279fbb1e5e4b2844ec9070808e34f0a8b834953201eab189bf6a929a4fffa3ebe062f072d04aad91b43c50fe5e29d95bf09b0ba8b8da52f4b02ca8302230865af2c2e1ce234db9fb25b82858a36d79aa43bf7f25772d9cde04365727118808f9183d447060d21fbdd0f7bae17cac2a6af02aea0e8a9c6668ab13b1ce4a03a2ec37cea050e36a7e58664640b12a1a6ee9ebb88bdc90e747334bcf27c1872692c22c81bd1785f1ae90d3e458f9e912322d9008ac9d689e23bfada3aa6ef6209313b9a7cab38790d9c2bd8eab2c813cc4e118935ace4cea31697cd956a5bc04f7fd86346a"},
    {13, "c4a8f7ff07b6ea6b1e99419766cd5d847fad53204fe8946152c1624873a6f65c096700c14402d26d444861b866c9759f5d4f05593aa9e95b9f6097eea55f4532c3ae16e06a262c88aa14890b5540ad0e8cff14de784a7f9ea450790371c4a36e3bcf604cb348d930f786c374128b39296ac362a5fd0d59c33059f369645cc3dd2fc2b5dad561420395f2e78c88e8874f2ad8cf8b1cfea9e66cd6e91c3ab49e29a9dd91339eba983cd7ca2ce9c2601e26dcac766c6f94773bf6e7f7d45024bd6603e96f8a3b06f333bc197fb29beeba2e8a085be33d4ddf2658b9f91cdd8959597eb8063af9bacc101a048e63c97950fa6b9707e782bab5e444091e654fda8918"},
    {14, "b85cb9833705ea853ab317013d17d28ee18211a438aa98d20a78830a615bde004d4bbe9d879aa7c4fe88bf79d8eec0440bf09341f6f9a502897f044fe1afe36ce96a43c588e555d66a4d243b5392ffb801e2a067fc5f06bac05d03fabb0a1accb6ebd6e585febb9dd88c569441d971dabcfd8bf10bdd02afb690f090a58a323410b1d879dbbd2526da091d80f63d0e9d13e98f9951e44f58b195908f7366fef564b229ef0dab7b33e0b3cd5b5bc347ffe37141d71aa8daa604c3073ef7838df73de97e53913cb99fb3f889b5e46eb78de3b06a166b84e1ee922bf5184498d6e067ce76b58aaa58ea3fe4fbb03e71a18525f174108cffea65b4aa6da85e396ea6"},
    {15, "97c0184d3181079c3e6da56f7d052e49d8da98150adde5719194409f2f95b6ac773419470598f5b23701b6ef76c93a6bd41764d3efe100ea3b2673287c331040d2577689c5f08180c821d7f599c2ae045b98f1e914c919f4816d289d461823baa4f7dbd9ca7a0ba3e8459785e1c088acaf2350e895def2222f0f6c6add79fe7c52c7b0734a61b4bb65d5f6682aaa0fb6ee28a4813d6ca60feaf3d9feb07dadc2d361f87598d8024b3a6ab11774cb0f4ba0a6d62c4576aca7af7f72686521f511a63110ef8ef05ca6b1b4f198dabd66335e22a5c70c34dd2ea681ec341f86685fb76c8d6bcf14524f0f8034de85affe323e29f8bbc1a21e71cd7475e5626a77c9"},
    {16, "9406abd0481c5ee592dafcf2b77cabb8e1667b5ae05293366b3b3e0ac07c83673fe13a378cf237540f54b7afb7a02bfbffb5424a0ba253fb86f25c64629d31742c4ac62987667a2cd6870cfcd404dc9b77b12a7eb2b7b6bb80f62de2868066ac9ad2071a6fb34ed56638b5c98f6189961d8fde7874d2e339fabcca96673355dc8decabec0e6f46526954f68f119b09691cab82166ebf043d4f4beb8a30da7caa5cb8ab86f66cc84bf0d0240ba4563d377e428d494226e8354f787ed2b4abacfbc877ad4586926ed67c6522968a1305d1b0e54b25a18368488b3b6e4a905cea5613722176db5d7ecb5b66596598bd5f75c455b22f16c7f010267320165f31916b"},
    {17, "ac4c90add54e1b8aaa542cbf84c5981ea2365ff011fb4f6b8b6b87d6931672e4bd4873e5e1290de03d21451ef7a8b0286c2c602d2d852a180fbd34ef502abe5ee8cb9bf2104aa084baa2f7b4a472c6c6b0adf3c36c242be16aa47c3b7014338483ed5a5753fc120b3ff86b17d04050a89b001376197e04bf05aca5a391f89acac4a89cdd06d5062e4dfab34c572e86e78facf10b812c2f6a67f4f998b2591a60da32c0daccc6c1ca195bb2feabefc3f3a6d7799f9b62182ed9c4aec20f0fd743b1ab71ebd9b00034a7feec208cdab095b6628a065cfbc8aab452bc6099231966f983ec734a628774b3e38613dab659199acf0f454786dc351de0ac305e502d15"},
    {18, "189dc6d3d27f93fefe32ba60a80e471091a0f22e5f6dc6f44d91ae71aaa4788f2838081d2dfbc83f73f26e684b27fecf87c24ff5ae3484822e1a14126647f9f3bbcb87c954004e7c27eddc3dc0745b9486643621ff694d340f2478128aa7c1a57245d899d3edbc85c60aa5479feade07e68c8272e1a6601bc67e17a8f51efd33e9e5ddcb7bcd1adb9e134d945ada563d5f12a6517ec2aa235e0b4da65c24fb7932e7aeb3e7614af3083ee787749d1d003662916b719594db669b40a14b247cbf59849dab19d930833a10e70e8096098fc8ce9f398734c2f6362179c518a66f63d0e1ed97bfca36d2644e28039639466036db172213c3885be103c7dfb3385a80"},
    {19, "44b8800ec43750fba19c7186af2653dd652af8d0bb72e6852a38c090afee801fed3aaefff62ae8d2c88a1a2afa9d5ea208d2f44ef9a93a20254fbbe2725fcbf5ac35c952149c23447c2aa81862ec58975d2679106326fdefc4d032549286e1d23de60822f8aa887aacc8719de38f4a64dec41d49cc01921e8cc415bf6315d9d3b0d57fdf861e03318617573e019de7e65a77ba9b4588bd8f7247f99801279926d88a04e2e260be49bc3e65a1189c126a376e58d0db437b154451f5b92748358a7533581e1e76f5158dc660833817c181b3b491a1e03481262ee178c3317e8b23bc8e3dd87b03a87857b1ce4d4a1fa469d7cd3d84505a5d8762ddfc522ac73418"},
    {20, "4ed421fa52dd1b828839b74b07f8ce6de0a16f752c893cf49dc55c23fd9ff3da5538df0979f42bc5b9a2dbb14c25b2da5a57ef814995ffba82441aac9633d537bddabd5f6229aff0cf92059aa3d22837e9daa06182968cc329f5926d53ac701d75add55494ced38da658faa7737aa19e49188b7a2621174e8897665dbe0d79c6f45b917112b5ff47b75aa45b4a8b1c0372185190e229e549fc58adfb1a452bf146bcbef3970737abd2e6a928baeb6c1b45f5f0a4188dc81e8d61f5f6413f3f5d3f03337eb8cf850e6388fd3a61fd33a831d552f6cecde6ef16f184dd91a0c717cfd81aa604277e799e46b09d0e55224503d7c9711596c88e7af73500afa30cf1"},
};

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
    const char *p;
    SQR_ERROR err;

    print_cpuid_info();

    printf("OpenSSL version: %s\n", OpenSSL_version(OPENSSL_VERSION));
    if( OpenSSL_version_num() < 0x10100000 ) {  // < 1.1.x ?
        printf("WARNING: OpenSSL version is too low (%08lx)\n", OpenSSL_version_num());
    }

    putchar('\n');

    err = sqr_allocate_state(RSA_2048_N, sizeof(RSA_2048_N), &state);
    if( err != SQR_ERR_OK )
        return 1;

    // test the lowest number this code supports
    err = sqr_calculate(state, 1, mhex, out, sizeof(out));
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000020e494579fb0b2420552a4b3140d5c1157f1d55c9506d9dfe161228dcf7d8e4b0c49a3ba23e0c58ff81163f977b09e5286244487c6c682e2324ffe821dd59000";
    if( strcmp(out, p) != 0 ) {
        printf("Failed in x^2:\n%s != \n%s\n", out, p);
        sqr_free_state(state);
        return 1;
    }

    err = sqr_calculate(state, 2, mhex, out, sizeof(out));
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000439f12e79f60df042c09accea1dea5dfd5031a690b773d31e49b8bce0109673a2d1216b9ceae9c909ca5e5903f51a410fc2d43385cae66b489950da94902063360f79e51b200729da14f9bc227785d450488f0b38edeeedb3b046bbcddaafab052b08ebc8e917a28c672054b49bd84a80c73d8e384e8b64fa2054c8f1000000";
    if( strcmp(out, p) != 0 ) {
        printf("Failed in x^2^2:\n%s != \n%s\n", out, p);
        sqr_free_state(state);
        return 1;
    }

    // this is the largest value in Montgomery for for our N: N-1 in Montgomery form
    err = sqr_calculate(state, 2, "a680e2f704919dc027ac89c4943ddcee74aae4c6c4b7eebfe070091df53500f8499288d4d25ec675be0b4b2f9c541f0e7dc61a605ef7616f1a5f83d7c4c85438cbcac933689bccd91dfd368b57b35ab9a75fca0cf1820d252ca941da4f8fd6a628f8198efa9c716e6535480af94c16d5d1b3a7d9dd29a8ce22ee31ea3615cb1e1b96290b6bf242b528fcc28c8f9ce9e5ab7a4d483dc14402245a58a02f999950c127c21c32a984bdaa615cef528708edeeb5e6cde0996b76f87741b6955594e6963314821d5af1b23e3f5a73f8b5094b368adfdc97eef388093d793193cbbdfc8dd2c191339e9238e6385dc7f8c1196aec42fbefed73fe9734f3fdce2a5a9fbb",
                        out, sizeof(out));
    p = "5411fae465c5275b443bb8ac114c055f99120e3792047c1673be69cc97ed107ca4906a375680935b328aa08f8ffcce2fe84710c4b537e6b7dd827d28d66dac18c235f828cbdbf161d727fc28682d99d7ee1aa83c00d2c1c2304aaa3ab7e63d9d5598d923fe208e8357858dbce1bf2f6f42f10548e9db7cffbe841e05edbaf44e6af43c800f9c1e91bff2b001e069c42db0fa70998e4bd427fe1ce33b685c1a650888c3d32900418e1c9b7c92b91b41381943bb2ed0e03ae29b1e7acd29146dbc521f16460c03a75586eb2210ed93bb21a73af0d9497b83a816d7f177b34ed66b85b887b9c9a646010b9e4fa645b94adb7a95f77d7a9518bf1c1b80974c8f4afe";
    if( strcmp(out, p) != 0 ) {
        printf("Failed in x^2^2:\n%s != \n%s\n", out, p);
        sqr_free_state(state);
        return 1;
    }

    // input x is close to sqrt(N): x^2 < N
    err = sqr_calculate(state, 2, "e20c000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", out, sizeof(out));
    p = "30848297497090cc5cd54f682fbf4993d295875dfbb4dfcde8caa1a9a1738b7e4f65ed8e02a840b8843522fe4066839cca8f10a04d96e137e1fec7d18e0d4687f20bff9e4568a005d60c9bb3f987e78142af0bc3835168ac44a47dd2d1c50273d6625ff734b05429d297abad01d8ef0ba33635d24703453cfe662f4ecd5996f3c3d57d654ddeb84710f63867684dbfabb31f2b5ed2293836d17bbcaab4d760ef52f963cffdc74116112b035ab67df1de8b1c608170326dc69e887664b2823c9362e8b1020cbfbfa0b97a9f71ccbab6118ecf376ca4800d37d8a7730a0850bb82331e2cd00897faeb87256995a965d4491d23bc8623247360ce8b1d11e29f9434";
    if( strcmp(out, p) != 0 ) {
        printf("Failed in x^2^2:\n%s != \n%s\n", out, p);
        sqr_free_state(state);
        return 1;
    }

    // test a reasonably large number
    err = sqr_calculate(state, 40, mhex, out, sizeof(out));
    p = "9d0073ae84cd3517d3db2f2b0000a82737248b576d8df421549d19320e7823a2c706dd9bdf21c4bf940aa870fd10465b3ea1f6131d629a4a849af14ee2afce19d905379963f93d39b19768ad23445a9a6d23111fad6335c1db3f5601f2151894e3034cac5f7c6f335f7a8a425805a27a152196f4437bc5424139a0fbd22201531e32adb71db66f5b9125e127cae4d9790a2174e81f1a019f6bf21fd9430a308f6bdd0241604f2bed98a37e4a169e15b376e3eb3f0e48cefaaad1fa17d63abdcdc173b72047a65652fe35ec421200ad72fcb13f492976cf93052ed741a14a2c343974e4dd1fb600329417fbd70e49139c562dfea58af68fdd1ab399639404d1ec";
    if( strcmp(out, p) != 0 ) {
        printf("Failed in x^2^40:\n%s != \n%s\n", out, p);
        sqr_free_state(state);
        return 1;
    }

    // even larger numbers from the kat array
    for(unsigned i = 0; i < sizeof(kat) / sizeof(kat[0]); i++) {
        err = sqr_calculate(state, 1 << kat[i].t, mhex, out, sizeof(out));
        if( strcmp(out, kat[i].y) != 0 ) {
            printf("Failed in x^2^2^%d:\n%s != \n%s\n", kat[i].t, out, kat[i].y);
            sqr_free_state(state);
            return 1;
        }
    }

    // finally a more relistic test with T = 2^t = 2^20
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

#ifndef SQR_H
#define SQR_H 1

typedef enum SQR_ERROR {
    SQR_ERR_OK = 0,
    SQR_ERR_OUT_OF_MEM = 1,
    SQR_ERR_BAD_PARAMS = 2,
    SQR_ERR_INTERNAL = 3
} SQR_ERROR;

typedef struct Sqr_State Sqr_State;

SQR_ERROR sqr_allocate_state( const uint8_t *N, unsigned Nl, Sqr_State **state );
SQR_ERROR sqr_free_state( Sqr_State *state );
SQR_ERROR sqr_calculate( Sqr_State *state, uint64_t T, const char *xhex, char *outhex, unsigned outl );

#endif // SQR_H

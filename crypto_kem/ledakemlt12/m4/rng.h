#ifndef RNG_H
#define RNG_H

#include <stddef.h>
#include <stdint.h>

#define RNG_SUCCESS     ( 0)
#define RNG_BAD_MAXLEN  (-1)
#define RNG_BAD_OUTBUF  (-2)
#define RNG_BAD_REQ_LEN (-3)
#define RNG_MAXLEN      (10 * 1024 * 1024)

typedef struct {
    uint8_t   buffer[16];
    size_t    buffer_pos;
    size_t    length_remaining;
    uint8_t   key[32];
    uint8_t   ctr[16];
} AES_XOF_struct;

int seedexpander(AES_XOF_struct *ctx, uint8_t *x, size_t xlen);
void seedexpander_from_trng(AES_XOF_struct *ctx, const uint8_t *trng_entropy);

#endif

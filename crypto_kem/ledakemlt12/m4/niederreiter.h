#ifndef NIEDERREITER_H
#define NIEDERREITER_H

#include "ring.h"
#include "params.h"
#include "rng.h"

typedef struct {
    uint8_t prng_seed[TRNG_BYTE_LENGTH];
    uint8_t rejections;
    uint8_t secondIterThreshold;
    uint8_t decryption_failure_secret[TRNG_BYTE_LENGTH];
} privateKeyNiederreiter_t;

typedef struct {
    DIGIT Mtr[(N0 - 1) * NUM_DIGITS_GF2X_ELEMENT];
} publicKeyNiederreiter_t;


void niederreiter_keygen(publicKeyNiederreiter_t *pk, privateKeyNiederreiter_t *sk);
void niederreiter_encrypt(DIGIT syndrome[], const publicKeyNiederreiter_t *pk, const DIGIT err[]);
int niederreiter_decrypt(DIGIT *err, const privateKeyNiederreiter_t *sk, const DIGIT *syndrome);

#endif

#include "api.h"
#include "niederreiter.h"
#include "randombytes.h"
#include "rng.h"
#include "utils.h"
#include "params.h"

#include <string.h>

/* Hacky workaround to achieve the same NIST KATs incase 
 * a different word size than w = 64 is used. Due to the 
 * representation, we are required to print the second digit
 * before the first to aquire the same ordering.
 * LEDAkemLT32, LEDAkemLT52, might require something different,
 * because there the number of digits is not even */
static void swap_digits(DIGIT *poly) {
    DIGIT t;
    for (size_t j = 0; j < POLY_DIGITS - 1; j += 2) {
        t = poly[j + 1];
        poly[j + 1] = poly[j];
        poly[j] = t;
    }
}

static void pack_error_inplace(DIGIT *error_digits) {
    for (size_t i = 0; i < N0; i++) {
        swap_digits(error_digits + i * POLY_DIGITS);
    }
}

/* IND-CCA2 Keygen */
int crypto_kem_keypair(uint8_t *pk, uint8_t *sk) {

    niederreiter_keygen((publicKeyNiederreiter_t *)pk, (privateKeyNiederreiter_t *) sk);
    swap_digits((DIGIT *)pk);

    return 0;
}

/* IND-CCA2 Encapsulation */
int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk) {
    AES_XOF_struct hashedAndTruncatedSeed_expander;
    POSITION_T errorPos[NUM_ERRORS_T];
    DIGIT error_vector[ERROR_DIGITS];
    uint8_t seed[TRNG_BYTE_LENGTH];
    uint8_t ss_input[2 * TRNG_BYTE_LENGTH] = {0};
    uint8_t hashedSeed[HASH_BYTE_LENGTH];
    uint8_t hashedAndTruncatedSeed[TRNG_BYTE_LENGTH] = {0};
    uint8_t hashedErrorVector[HASH_BYTE_LENGTH];
    uint8_t hashedAndTruncatedErrorVector[TRNG_BYTE_LENGTH] = {0};
    uint8_t maskedSeed[TRNG_BYTE_LENGTH];

    swap_digits((DIGIT *) pk);

    randombytes(seed, TRNG_BYTE_LENGTH);
    memcpy(ss_input, seed, TRNG_BYTE_LENGTH);

    HASH_FUNCTION(ss, ss_input, 2 * TRNG_BYTE_LENGTH);
    HASH_FUNCTION(hashedSeed, seed, TRNG_BYTE_LENGTH);

    memcpy(hashedAndTruncatedSeed, hashedSeed, TRNG_BYTE_LENGTH);

    memset(&hashedAndTruncatedSeed_expander, 0x00, sizeof(AES_XOF_struct));
    seedexpander_from_trng(&hashedAndTruncatedSeed_expander, hashedAndTruncatedSeed);
    rand_error_pos(errorPos, &hashedAndTruncatedSeed_expander);
    expand_error(error_vector, errorPos);

    niederreiter_encrypt((DIGIT *) ct, (const publicKeyNiederreiter_t *) pk, error_vector);

    pack_error_inplace(error_vector);
    HASH_FUNCTION(hashedErrorVector, (const uint8_t *)error_vector, ERROR_BYTES);

    memcpy(hashedAndTruncatedErrorVector, hashedErrorVector, TRNG_BYTE_LENGTH);

    for (size_t i = 0; i < TRNG_BYTE_LENGTH; ++i) {
        maskedSeed[i] = seed[i] ^ hashedAndTruncatedErrorVector[i];
    }
    
    swap_digits((DIGIT *) pk);
    swap_digits((DIGIT *) ct);
    memcpy(ct + (POLY_BYTES), maskedSeed, TRNG_BYTE_LENGTH);

    return 0;
}


/* IND-CCA2 Decapsulation  */
int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk) {
    AES_XOF_struct hashedAndTruncatedSeed_expander;
    POSITION_T reconstructed_errorPos[NUM_ERRORS_T];
    DIGIT reconstructed_error_vector[ERROR_DIGITS];
    DIGIT decoded_error_vector[ERROR_DIGITS];
    uint8_t hashedErrorVector[HASH_BYTE_LENGTH];
    uint8_t hashedAndTruncatedErrorVector[TRNG_BYTE_LENGTH] = {0};
    uint8_t decoded_seed[TRNG_BYTE_LENGTH];
    uint8_t hashed_decoded_seed[HASH_BYTE_LENGTH];
    uint8_t hashedAndTruncated_decoded_seed[TRNG_BYTE_LENGTH] = {0};
    uint8_t ss_input[2 * TRNG_BYTE_LENGTH], tail[TRNG_BYTE_LENGTH] = {0};
    int decode_ok, decrypt_ok, equal;

    swap_digits((DIGIT *) ct);

    decode_ok = niederreiter_decrypt(decoded_error_vector,
                                     (const privateKeyNiederreiter_t *)sk,
                                     (const DIGIT *) ct);

    pack_error_inplace(decoded_error_vector);
    HASH_FUNCTION(hashedErrorVector, (const uint8_t *) decoded_error_vector, ERROR_BYTES);

    memcpy(hashedAndTruncatedErrorVector, hashedErrorVector, TRNG_BYTE_LENGTH);

    for (size_t i = 0; i < TRNG_BYTE_LENGTH; ++i) {
        decoded_seed[i] = ct[(POLY_BYTES) + i] ^ hashedAndTruncatedErrorVector[i];
    }

    HASH_FUNCTION(hashed_decoded_seed, decoded_seed, TRNG_BYTE_LENGTH);

    memcpy(hashedAndTruncated_decoded_seed, hashed_decoded_seed, TRNG_BYTE_LENGTH);

    memset(&hashedAndTruncatedSeed_expander, 0x00, sizeof(AES_XOF_struct));
    seedexpander_from_trng(&hashedAndTruncatedSeed_expander, hashed_decoded_seed);

    rand_error_pos(reconstructed_errorPos, &hashedAndTruncatedSeed_expander);

    expand_error(reconstructed_error_vector, reconstructed_errorPos);

    pack_error_inplace(reconstructed_error_vector);
    equal = gf2x_verify(decoded_error_vector, reconstructed_error_vector, ERROR_DIGITS);
    // equal == 0, if the reconstructed error vector match !!!

    decrypt_ok = (decode_ok == 1 && equal == 0);

    memcpy(ss_input, decoded_seed, TRNG_BYTE_LENGTH);
    memcpy(ss_input + sizeof(decoded_seed), tail, TRNG_BYTE_LENGTH);

    // Overwrite on failure
    cmov(ss_input + sizeof(decoded_seed),
         ((const privateKeyNiederreiter_t *) sk)->decryption_failure_secret,
         TRNG_BYTE_LENGTH, !decrypt_ok);

    HASH_FUNCTION(ss, ss_input, 2 * TRNG_BYTE_LENGTH);

    swap_digits((DIGIT *) ct);

    return 0;
}

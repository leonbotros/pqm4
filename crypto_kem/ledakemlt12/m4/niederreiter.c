#include "hq_gen.h"
#include "bf_decoding.h"
#include "dfr_test.h"
#include "ring.h"
#include "niederreiter.h"
#include "params.h"
#include "randombytes.h"
#include "rng.h"
#include "utils.h"

#include <string.h>

void niederreiter_keygen(publicKeyNiederreiter_t *pk, privateKeyNiederreiter_t *sk) {
    AES_XOF_struct keys_expander;
    POSITION_T HPosOnes[N0][DV];
    POSITION_T QPosOnes[N0][M];
    POSITION_T LPosOnes[N0][DV * M];
    POSITION_T auxPosOnes[DV * M];
    uint8_t processedQOnes[N0];
    DIGIT* Ml = pk->Mtr; // RAM reduction, does not work for n0 != 2
    DIGIT* Ln0dense = Ml;
    DIGIT* Ln0Inv = Ln0dense;
    int is_L_full;
    int isDFRok = 0;

    memset(&keys_expander, 0x00, sizeof(AES_XOF_struct));
    randombytes(sk->prng_seed, TRNG_BYTE_LENGTH);
    seedexpander_from_trng(&keys_expander, sk->prng_seed);

    sk->rejections = (uint8_t) 0;
    do {
        generateHPosOnes(HPosOnes, &keys_expander);
        generateQPosOnes(QPosOnes, &keys_expander);
        for (int i = 0; i < N0; i++) {
            for (int j = 0; j < DV * M; j++) {
                LPosOnes[i][j] = INVALID_POS_VALUE;
            }
        }

        memset(processedQOnes, 0x00, sizeof(processedQOnes));
        for (int colQ = 0; colQ < N0; colQ++) {
            for (int i = 0; i < N0; i++) {
                gf2x_mod_mul_sparse(DV * M, auxPosOnes,
                                    DV, HPosOnes[i],
                                    qBlockWeights[i][colQ], QPosOnes[i] + processedQOnes[i]);
                gf2x_mod_add_sparse(DV * M, LPosOnes[colQ],
                                    DV * M, LPosOnes[colQ],
                                    DV * M, auxPosOnes);
                processedQOnes[i] += qBlockWeights[i][colQ];
            }
        }
        is_L_full = 1;
        for (size_t i = 0; i < N0; i++) {
            is_L_full = is_L_full && (LPosOnes[i][DV * M - 1] != INVALID_POS_VALUE);
        }
        sk->rejections = sk->rejections + 1;
        if (is_L_full) {
            isDFRok = DFR_test(LPosOnes, &(sk->secondIterThreshold));
        }
    } while (!is_L_full || !isDFRok);
    sk->rejections = sk->rejections - 1;

    seedexpander(&keys_expander, sk->decryption_failure_secret, TRNG_BYTE_LENGTH);

    memset(Ln0dense, 0x00, POLY_BYTES);
    for (size_t j = 0; j < DV * M; j++) {
        if (LPosOnes[N0 - 1][j] != INVALID_POS_VALUE) {
            gf2x_set_coeff(Ln0dense, LPosOnes[N0 - 1][j], 1);
        }
    }

    gf2x_mod_inverse(Ln0Inv, Ln0dense);
    for (size_t i = 0; i < N0 - 1; i++) {
        gf2x_mod_mul_dense_to_sparse(Ml, Ln0Inv, LPosOnes[i], DV * M);
    }

    for (size_t i = 0; i < N0 - 1; i++) {
        gf2x_transpose_in_place(Ml);
    }
}

void niederreiter_encrypt(DIGIT syndrome[],
                          const publicKeyNiederreiter_t *pk,
                          const DIGIT err[]) {
    DIGIT tmp[POLY_DIGITS];

    memset(syndrome, 0x00, POLY_BYTES);
    for (size_t i = 0; i < N0 - 1; i++) {
        gf2x_mod_mul(tmp, pk->Mtr + i * POLY_DIGITS, err + i * POLY_DIGITS);
        gf2x_mod_add(syndrome, syndrome, tmp);
    }
    gf2x_mod_add(syndrome, syndrome, err + (N0 - 1) * POLY_DIGITS);
}

int niederreiter_decrypt(DIGIT *err, const privateKeyNiederreiter_t *sk, const DIGIT *syndrome) {
    AES_XOF_struct niederreiter_decrypt_expander;
    POSITION_T HPosOnes[N0][DV];
    POSITION_T QPosOnes[N0][M];
    POSITION_T LPosOnes[N0][DV * M];
    POSITION_T auxPosOnes[DV * M];
    POSITION_T HtrPosOnes[N0][DV];
    POSITION_T QtrPosOnes[N0][M];
    POSITION_T auxSparse[DV * M];
    POSITION_T Ln0trSparse[DV * M];
    DIGIT privateSyndrome[POLY_DIGITS];
    uint8_t processedQOnes[N0];
    uint8_t mask;
    int rejections = sk->rejections;
    int decrypt_ok = 0;
    int err_weight;

    seedexpander_from_trng(&niederreiter_decrypt_expander, sk->prng_seed);
    do {
        generateHPosOnes(HPosOnes, &niederreiter_decrypt_expander);
        generateQPosOnes(QPosOnes, &niederreiter_decrypt_expander);

        for (size_t i = 0; i < N0; i++) {
            for (size_t j = 0; j < DV * M; j++) {
                LPosOnes[i][j] = INVALID_POS_VALUE;
            }
        }

        memset(processedQOnes, 0x00, sizeof(processedQOnes));
        for (size_t colQ = 0; colQ < N0; colQ++) {
            for (size_t i = 0; i < N0; i++) {
                gf2x_mod_mul_sparse(DV * M, auxPosOnes,
                                    DV, HPosOnes[i],
                                    qBlockWeights[i][colQ], QPosOnes[i] + processedQOnes[i]);
                gf2x_mod_add_sparse(DV * M, LPosOnes[colQ],
                                    DV * M, LPosOnes[colQ],
                                    DV * M, auxPosOnes);
                processedQOnes[i] += qBlockWeights[i][colQ];
            }
        }
        rejections--;
    } while (rejections >= 0);

    transposeHPosOnes(HtrPosOnes, HPosOnes);
    transposeQPosOnes(QtrPosOnes, QPosOnes);

    for (size_t i = 0; i < DV * M; i++) {
        Ln0trSparse[i] = INVALID_POS_VALUE;
        auxSparse[i] = INVALID_POS_VALUE;
    }

    for (size_t i = 0; i < N0; i++) {
        gf2x_mod_mul_sparse(DV * M, auxSparse,
                            DV, HPosOnes[i],
                            qBlockWeights[i][N0 - 1], &QPosOnes[i][M - qBlockWeights[i][N0 - 1]]);
        gf2x_mod_add_sparse(DV * M, Ln0trSparse,
                            DV * M, Ln0trSparse,
                            DV * M, auxSparse);
    }
    gf2x_transpose_in_place_sparse(DV * M, Ln0trSparse);

    gf2x_mod_mul_dense_to_sparse(privateSyndrome, syndrome, Ln0trSparse, DV * M);

    
    decrypt_ok = bf_decoding(err,
                             (const POSITION_T (*)[DV]) HtrPosOnes,
                             (const POSITION_T (*)[M])  QtrPosOnes,
                             privateSyndrome, sk->secondIterThreshold);

    err_weight = 0;
    for (size_t i = 0 ; i < N0; i++) {
        err_weight += population_count(err + (POLY_DIGITS * i));
    }
    decrypt_ok = decrypt_ok && (err_weight == NUM_ERRORS_T);

    /* overwrite err with mockup on decryption failure, saves ram */
    gf2x_cmov(err, syndrome, POLY_DIGITS, !decrypt_ok);
    cmov((uint8_t *)err + POLY_BYTES, sk->decryption_failure_secret, TRNG_BYTE_LENGTH, !decrypt_ok);
    
    // conditionally pad with remainder with zeroes
    mask = -(uint8_t) decrypt_ok;
    uint8_t *err_bytes = (uint8_t *) err;
    for(size_t i = 0; i < (N0 - 1) * POLY_BYTES - TRNG_BYTE_LENGTH; i++) {
        err_bytes[POLY_BYTES + TRNG_BYTE_LENGTH + i] = (err_bytes[POLY_BYTES + TRNG_BYTE_LENGTH + i] & mask) ^ (0x00 & ~mask);
    }

    return decrypt_ok;
}

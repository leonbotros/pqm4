#ifndef GF2X_ARITH_MOD_XPLUSONE_H
#define GF2X_ARITH_MOD_XPLUSONE_H

#include "params.h"

#include "gf2x.h"
#include "rng.h"

void gf2x_copy(DIGIT dest[], const DIGIT in[], size_t n);
DIGIT gf2x_get_coeff(const DIGIT poly[], size_t exponent);
void gf2x_set_coeff(DIGIT poly[], size_t exponent, DIGIT value);
void gf2x_toggle_coeff(DIGIT poly[], size_t exponent);
int population_count(const DIGIT *poly);
void gf2x_mod_add(DIGIT Res[], const DIGIT A[], const DIGIT B[]);
void gf2x_mod(DIGIT out[],  const DIGIT in[]);
void gf2x_mod_mul(DIGIT Res[], const DIGIT A[], const DIGIT B[]);
void gf2x_transpose_in_place(DIGIT A[]);
void rand_circulant_sparse_block(POSITION_T *pos_ones, size_t countOnes, AES_XOF_struct *seed_expander_ctx);
void rand_circulant_blocks_sequence(DIGIT *sequence, AES_XOF_struct *seed_expander_ctx);
void rand_error_pos(POSITION_T errorPos[NUM_ERRORS_T], AES_XOF_struct *seed_expander_ctx);
void expand_error(DIGIT sequence[ERROR_DIGITS], const POSITION_T errorPos[NUM_ERRORS_T]);
void gf2x_mod_add_sparse(size_t sizeR, POSITION_T Res[], size_t sizeA, const POSITION_T A[], size_t sizeB, const POSITION_T B[]);
void gf2x_transpose_in_place_sparse(size_t sizeA, POSITION_T A[]);
int gf2x_mod_inverse(DIGIT out[], const DIGIT in[]);
void gf2x_mod_mul_sparse(size_t sizeR, POSITION_T Res[], size_t sizeA,  const POSITION_T A[], size_t sizeB, const POSITION_T B[]);
void gf2x_mod_mul_dense_to_sparse(DIGIT *R, const DIGIT *dense, const POSITION_T *sparse, size_t n);
void gf2x_tobytes(uint8_t *bytes, const DIGIT *poly);
void gf2x_frombytes(DIGIT *poly, const uint8_t *poly_bytes);

#endif

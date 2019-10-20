#ifndef GF2X_ARITH_H
#define GF2X_ARITH_H

#include "params.h"

void gf2x_add(DIGIT Res[], const DIGIT A[], const DIGIT B[], size_t n);
void gf2x_cmov(DIGIT *r, const DIGIT *a, size_t len, int c);
void right_bit_shift_n(size_t length, DIGIT in[], size_t amount);
void left_bit_shift_n(size_t length, DIGIT in[], size_t amount);
void gf2x_mul(DIGIT *R, const DIGIT *A, const DIGIT *B, size_t n);

#endif

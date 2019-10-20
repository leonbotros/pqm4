#ifndef MULL_SMALL_H
#define MULL_SMALL_H

#include "gf2x.h"

extern void asm_gf2x_mul1(DIGIT *, DIGIT, DIGIT);
#define gf2x_mul1 asm_gf2x_mul1

//void gf2x_mul1(DIGIT *c, const DIGIT a, const DIGIT b);
void gf2x_mul2(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul3(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul4(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul5(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul6(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul7(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul8(DIGIT *c, const DIGIT *a, const DIGIT *b);
void gf2x_mul9(DIGIT *c, const DIGIT *a, const DIGIT *b);

#endif

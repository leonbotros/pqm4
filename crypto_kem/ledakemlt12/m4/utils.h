#ifndef UTILS_H
#define UTILS_H

#include "gf2x.h"

#include <stdint.h>

int gf2x_verify(const DIGIT *a, const DIGIT *b, size_t len);
void cmov(uint8_t *r, const uint8_t *a, size_t len, int cond);

#endif

#include "ring.h"
#include "rng.h"
#include "sort.h"

#include <string.h>

/* assembly subroutines */
extern void asm_gf2x_shift_right(DIGIT *);
extern void asm_gf2x_rotate_right(DIGIT *);
extern void asm_gf2x_mult_acc_scalar(DIGIT *, DIGIT *, DIGIT);
extern void asm_gf2x_cswap(DIGIT *, DIGIT *, DIGIT);

void gf2x_copy(DIGIT dest[], const DIGIT in[], size_t n) {
    for (size_t i = 0; i < n; i++) {
        dest[i] = in[i];
    }
}

/* returns the coefficient of the x^exponent term as the LSB of a digit */
DIGIT gf2x_get_coeff(const DIGIT poly[], size_t exponent) {
    size_t straightIdx = (POLY_DIGITS * DIGIT_SIZE_b - 1) - exponent;
    size_t digitIdx = straightIdx / DIGIT_SIZE_b;
    size_t inDigitIdx = straightIdx % DIGIT_SIZE_b;
    return (poly[digitIdx] >> (DIGIT_SIZE_b - 1 - inDigitIdx)) & ((DIGIT) 1) ;
}

/* sets the coefficient of the x^exponent term as the LSB of a digit */
void gf2x_set_coeff(DIGIT poly[], size_t exponent, DIGIT value) {
    size_t straightIdx = (POLY_DIGITS * DIGIT_SIZE_b - 1) - exponent;
    size_t digitIdx = straightIdx / DIGIT_SIZE_b;
    size_t inDigitIdx = straightIdx % DIGIT_SIZE_b;

    /* clear given coefficient */
    DIGIT mask = ~(((DIGIT) 1) << (DIGIT_SIZE_b - 1 - inDigitIdx));
    poly[digitIdx] = poly[digitIdx] & mask;
    poly[digitIdx] = poly[digitIdx] | ((value & ((DIGIT) 1)) << (DIGIT_SIZE_b - 1 - inDigitIdx));
}

/* toggles (flips) the coefficient of the x^exponent term as the LSB of a digit */
void gf2x_toggle_coeff(DIGIT poly[], size_t exponent) {
    size_t straightIdx = (POLY_DIGITS * DIGIT_SIZE_b - 1) - exponent;
    size_t digitIdx = straightIdx / DIGIT_SIZE_b;
    size_t inDigitIdx = straightIdx % DIGIT_SIZE_b;

    /* clear given coefficient */
    DIGIT mask = (((DIGIT) 1) << (DIGIT_SIZE_b - 1 - inDigitIdx));
    poly[digitIdx] = poly[digitIdx] ^ mask;
}

/* population count for an unsigned 32-bit integer
   Source: Hacker's delight, p.66  */
static int popcount_uint32t(uint32_t x) {
    x -= (x >> 1) & 0x55555555;
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0f0f0f0f;
    return (int)(x * 0x01010101) >> 24;
}

/* population count for a single polynomial */
int population_count(const DIGIT *poly) {
    int ret = 0;
    for (int i = POLY_DIGITS - 1; i >= 0; i--) {
        ret += popcount_uint32t(poly[i]);
    }
    return ret;
}

void gf2x_mod_add(DIGIT Res[], const DIGIT A[], const DIGIT B[]) {
    gf2x_add(Res, A, B, POLY_DIGITS);
}

void gf2x_mod(DIGIT out[],  const DIGIT in[]) {
    DIGIT aux[POLY_DIGITS + 1];

    memcpy(aux, in, (POLY_DIGITS + 1) * DIGIT_SIZE_B);
    right_bit_shift_n(POLY_DIGITS + 1, aux, TAIL);
    gf2x_add(out, aux + 1, in + POLY_DIGITS, POLY_DIGITS);
    out[0] &= TAIL_MASK;
}

/* asm to reverse bits of digit */
static DIGIT reverse_digit(DIGIT x) {
    asm("rbit %[rd], %[rm]" : [rd] "=r" (x) : [rm] "r" (x));
    return x;
}

void gf2x_transpose_in_place(DIGIT A[]) {
    /* it keeps the lsb in the same position and
     * inverts the sequence of the remaining bits */

    DIGIT mask = (DIGIT)0x1;
    DIGIT rev1, rev2, a00;
    int slack_bits_amount = POLY_DIGITS * DIGIT_SIZE_b - P;

    a00 = A[POLY_DIGITS - 1] & mask;
    asm_gf2x_shift_right(A);

    for (size_t i = POLY_DIGITS - 1; i >= (POLY_DIGITS + 1) / 2; i--) {
        rev1 = reverse_digit(A[i]);
        rev2 = reverse_digit(A[POLY_DIGITS - 1 - i]);
        A[i] = rev2;
        A[POLY_DIGITS - 1 - i] = rev1;
    }

#if(POLY_DIGITS % 2 == 1)
   A[POLY_DIGITS / 2] = reverse_digit(A[POLY_DIGITS / 2]);
#endif

    if (slack_bits_amount) {
        right_bit_shift_n(POLY_DIGITS, A, slack_bits_amount);
    }
    A[POLY_DIGITS - 1] = (A[POLY_DIGITS - 1] & (~mask)) | a00;
}

/* returns -1 mask if x != 0, otherwise 0 */
static inline int32_t nonzero(DIGIT x) {
    DIGIT t = x;
    t = (~t) + 1;
    t >>= DIGIT_SIZE_b - 1;
    return -((int32_t)t);
}

/* returns -1 mask if x < 0 else 0 */
static inline int32_t negative(int x) {
    uint32_t u = x;
    u >>= 31;
    return -((int32_t)u);
}

/* return f(0) as digit */
static inline DIGIT lsb(const DIGIT *p) {
    DIGIT mask = (DIGIT)1;
    return p[POLY_DIGITS - 1] & mask;
}

/* constant-time inverse, source: gcd.cr.yp.to */
int gf2x_mod_inverse(DIGIT out[], const DIGIT in[]) {
    int32_t swap, delta = 1;
    DIGIT g0_mask;

    DIGIT f[NUM_DIGITS_GF2X_MODULUS] = {0}; // f = x^P + 1
    DIGIT g[POLY_DIGITS];       // g = in
    DIGIT *v = out;                         // v = 0, save space
    DIGIT r[POLY_DIGITS] = {0}; // r = 1

    f[NUM_DIGITS_GF2X_MODULUS - 1] = 1;
    f[0] |= ((DIGIT)1 << MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS);

    for (size_t i = 0; i < POLY_DIGITS; i++) {
        g[i] = in[i];
    }

    for (size_t i = 0; i < POLY_DIGITS; i++) {
        v[i] = 0;
    }

    r[POLY_DIGITS - 1] = 1;

    for (int loop = 0; loop < 2 * P - 1; ++loop) {

        swap = negative(-delta) & nonzero(lsb(g));              // swap = -1 if -delta < 0 AND g(0) != 0
        delta ^= swap & (delta ^ -delta);                       // cond swap delta with -delta if swap
        delta++;

        asm_gf2x_cswap(f, g, swap);
        asm_gf2x_cswap(v, r, swap);

        g0_mask = ~lsb(g) + 1;

        // g = (g - g0 * f) / x
        asm_gf2x_mult_acc_scalar(g, f, g0_mask);
        asm_gf2x_shift_right(g);

        // r = (r - g0 * v) / x
        asm_gf2x_mult_acc_scalar(r, v, g0_mask);
        asm_gf2x_rotate_right(r);

    }

    return nonzero(delta) + 1; // 0 if fail, 1 if success
}

void gf2x_mod_mul(DIGIT Res[], const DIGIT A[], const DIGIT B[]) {
    DIGIT aux[2 * POLY_DIGITS];
    
    gf2x_mul(aux, A, B, POLY_DIGITS);
    gf2x_mod(Res, aux);
}

/* constant-time sparse-dense multiplication 
 * assumes no INVALID_POS in sparse, does not require
 * sorted sparse representaton */
void gf2x_mod_mul_dense_to_sparse(DIGIT *R, const DIGIT *dense,
                                   const POSITION_T *sparse, size_t n) {
    DIGIT tmpR[RLEN]; 
    DIGIT p0[RLEN];
    DIGIT p1[RLEN];
    DIGIT *v = p0;
    DIGIT *w = p1;
    DIGIT *ptr;
    size_t s, ds, ids;  // total shift, digit_shift and indigit_shift
    size_t i, j;
    int8_t b;           // bit position in pos in binary representation of index
    DIGIT mask;

    for (i = 0; i < n; i++) {
        s = sparse[i];
        
        for (j = 0; j < POLY_DIGITS; j++) {
            w[j] = 0;
            w[POLY_DIGITS + j] = dense[j];
        }

        for (b = P_BITS - 1; b >= LOG_DIGIT_SIZE_b; b--) {
            ptr = v; v = w; w = ptr; // use output of last loop as input of next loop

            mask = (DIGIT) -((s >> b) & 1);     // all-zero or all-one depending on bit 
            ds = 1 << (b - LOG_DIGIT_SIZE_b);   // 2^b / 2^5 = 2^(b-5) number of digits to rotate 
         
            for (j = 0; j < RLEN - ds; j++) {
                w[j] = (v[j + ds] & mask) ^ (v[j] & ~mask);
            }
            for(; j < RLEN; j++) {
                w[j] = (((v[j + ds - RLEN] << TAIL) | (v[j + ds - RLEN + 1] >> (DIGIT_SIZE_b - TAIL))) & mask) ^ (v[j] & ~mask);
            }
        }
       
        ids = s & ((1 << LOG_DIGIT_SIZE_b) - 1); // ids = 0 .. 31
        left_bit_shift_n(RLEN, w, ids); // shift inside digits by ids

        if (i == 0)
            gf2x_copy(tmpR, w, RLEN);
        else
            gf2x_add(tmpR, tmpR, w, RLEN);
    }
    
    gf2x_mod(R, tmpR);
}

void gf2x_transpose_in_place_sparse(size_t sizeA, POSITION_T A[]) {
    POSITION_T t;
    size_t i = 0, j;

    if (A[i] == 0) {
        i = 1;
    }
    j = i;

    for (; i < sizeA && A[i] != INVALID_POS_VALUE; i++) {
        A[i] = P - A[i];
    }

    for (i -= 1; j < i; j++, i--) {
        t = A[j];
        A[j] = A[i];
        A[i] = t;
    }

}

void gf2x_mod_mul_sparse(size_t sizeR, POSITION_T Res[],
                         size_t sizeA, const POSITION_T A[],
                         size_t sizeB, const POSITION_T B[]) {

    POSITION_T prod;
    POSITION_T lastReadPos;
    size_t duplicateCount;
    size_t write_idx, read_idx;

    /* compute all the coefficients, filling invalid positions with P*/
    size_t lastFilledPos = 0;
    for (size_t i = 0 ; i < sizeA ; i++) {
        for (size_t j = 0 ; j < sizeB ; j++) {
            prod = A[i] + B[j];
            prod = ( (prod >= P) ? prod - P : prod);
            if ((A[i] != INVALID_POS_VALUE) &&
                    (B[j] != INVALID_POS_VALUE)) {
                Res[lastFilledPos] = prod;
            } else {
                Res[lastFilledPos] = INVALID_POS_VALUE;
            }
            lastFilledPos++;
        }
    }
    while (lastFilledPos < sizeR) {
        Res[lastFilledPos] = INVALID_POS_VALUE;
        lastFilledPos++;
    }

    uint32_sort(Res, sizeR);

    /* eliminate duplicates */
    write_idx = read_idx = 0;
    while (read_idx < sizeR  && Res[read_idx] != INVALID_POS_VALUE) {
        lastReadPos = Res[read_idx];
        read_idx++;
        duplicateCount = 1;
        while ( (Res[read_idx] == lastReadPos) && (Res[read_idx] != INVALID_POS_VALUE)) {
            read_idx++;
            duplicateCount++;
        }
        if (duplicateCount % 2) {
            Res[write_idx] = lastReadPos;
            write_idx++;
        }
    }
    /* fill remaining cells with INVALID_POS_VALUE */
    for (; write_idx < sizeR; write_idx++) {
        Res[write_idx] = INVALID_POS_VALUE;
    }
}

/* the implementation is safe even in case A or B alias with the result
 * PRE: A and B should be sorted, disjunct arrays ending with INVALID_POS_VALUE */
void gf2x_mod_add_sparse(size_t sizeR, POSITION_T Res[],
                         size_t sizeA, const POSITION_T A[],
                         size_t sizeB, const POSITION_T B[]) {

    POSITION_T tmpRes[DV * M];
    size_t idxA = 0, idxB = 0, idxR = 0;
    while ( idxA < sizeA  &&
            idxB < sizeB  &&
            A[idxA] != INVALID_POS_VALUE &&
            B[idxB] != INVALID_POS_VALUE ) {

        if (A[idxA] == B[idxB]) {
            idxA++;
            idxB++;
        } else {
            if (A[idxA] < B[idxB]) {
                tmpRes[idxR] = A[idxA];
                idxA++;
            } else {
                tmpRes[idxR] = B[idxB];
                idxB++;
            }
            idxR++;
        }
    }

    while (idxA < sizeA && A[idxA] != INVALID_POS_VALUE) {
        tmpRes[idxR] = A[idxA];
        idxA++;
        idxR++;
    }

    while (idxB < sizeB && B[idxB] != INVALID_POS_VALUE) {
        tmpRes[idxR] = B[idxB];
        idxB++;
        idxR++;
    }

    while (idxR < sizeR) {
        tmpRes[idxR] = INVALID_POS_VALUE;
        idxR++;
    }
    memcpy(Res, tmpRes, sizeof(POSITION_T)*sizeR);

}

/* Return a uniform random value in the range 0..n-1 inclusive,
 * applying a rejection sampling strategy and exploiting as a random source
 * the NIST seedexpander seeded with the proper key.
 * Assumes that the maximum value for the range n is 2^32-1
 */
static uint32_t rand_range(const unsigned int n, const int logn, AES_XOF_struct *seed_expander_ctx) {
    unsigned long required_rnd_bytes = (logn + 7) / 8;
    unsigned char rnd_char_buffer[4];
    uint32_t rnd_value;
    uint32_t mask = ( (uint32_t)1 << logn) - 1;

    do {
        seedexpander(seed_expander_ctx, rnd_char_buffer, required_rnd_bytes);
        /* obtain an endianness independent representation of the generated random
         bytes into an unsigned integer */
        rnd_value =  ((uint32_t)rnd_char_buffer[3] << 24) +
                     ((uint32_t)rnd_char_buffer[2] << 16) +
                     ((uint32_t)rnd_char_buffer[1] <<  8) +
                     ((uint32_t)rnd_char_buffer[0] <<  0) ;
        rnd_value = mask & rnd_value;
    } while (rnd_value >= n);

    return rnd_value;
}

/* Obtains fresh randomness and seed-expands it until all the required positions
 * for the '1's in the circulant block are obtained */
void rand_circulant_sparse_block(POSITION_T *pos_ones,
                                 size_t countOnes,
                                 AES_XOF_struct *seed_expander_ctx) {

    size_t duplicated, placedOnes = 0;
    POSITION_T p;

    while (placedOnes < countOnes) {
        p = rand_range(NUM_BITS_GF2X_ELEMENT,
                       P_BITS,
                       seed_expander_ctx);
        duplicated = 0;
        for (size_t j = 0; j < placedOnes; j++) {
            if (pos_ones[j] == p) {
                duplicated = 1;
            }
        }
        if (duplicated == 0) {
            pos_ones[placedOnes] = p;
            placedOnes++;
        }
    }
}

/* Returns random weight-t circulant block */
void rand_circulant_blocks_sequence(DIGIT sequence[ERROR_DIGITS],
                                    AES_XOF_struct *seed_expander_ctx) {

    size_t polyIndex, duplicated, counter = 0;
    POSITION_T p, exponent, rndPos[NUM_ERRORS_T];

    memset(sequence, 0x00, ERROR_BYTES);

    while (counter < NUM_ERRORS_T) {
        p = rand_range(N0 * NUM_BITS_GF2X_ELEMENT, P_BITS, seed_expander_ctx);
        duplicated = 0;
        for (size_t j = 0; j < counter; j++) {
            if (rndPos[j] == p) {
                duplicated = 1;
            }
        }
        if (duplicated == 0) {
            rndPos[counter] = p;
            counter++;
        }
    }
    for (size_t j = 0; j < counter; j++) {
        polyIndex = rndPos[j] / P;
        exponent = rndPos[j] % P;
        gf2x_set_coeff( sequence + POLY_DIGITS * polyIndex, exponent,
                        ( (DIGIT) 1));
    }

}


void rand_error_pos(POSITION_T errorPos[NUM_ERRORS_T],
                    AES_XOF_struct *seed_expander_ctx) {

    int duplicated;
    size_t counter = 0;

    while (counter < NUM_ERRORS_T) {
        POSITION_T p = rand_range(N0 * NUM_BITS_GF2X_ELEMENT, P_BITS, seed_expander_ctx);
        duplicated = 0;
        for (size_t j = 0; j < counter; j++) {
            if (errorPos[j] == p) {
                duplicated = 1;
            }
        }
        if (duplicated == 0) {
            errorPos[counter] = p;
            counter++;
        }
    }
}

void expand_error(DIGIT sequence[ERROR_DIGITS],
                  const POSITION_T errorPos[NUM_ERRORS_T]) {

    size_t polyIndex;
    POSITION_T exponent;

    memset(sequence, 0x00, ERROR_BYTES);
    for (int j = 0; j < NUM_ERRORS_T; j++) {
        polyIndex = errorPos[j] / P;
        exponent = errorPos[j] % P;
        gf2x_set_coeff( sequence + POLY_DIGITS * polyIndex, exponent,
                        ( (DIGIT) 1));
    }
}


void gf2x_tobytes(uint8_t *bytes, const DIGIT *poly) {
    size_t i, j;
    for (i = 0; i < POLY_DIGITS; i++) {
        for (j = 0; j < DIGIT_SIZE_B; j++) {
            bytes[i * DIGIT_SIZE_B + j] = (uint8_t) (poly[i] >> 8 * j);
        }
    }
}

void gf2x_frombytes(DIGIT *poly, const uint8_t *poly_bytes) {
    size_t i, j;
    for (i = 0; i < POLY_DIGITS; i++) {
        poly[i] = (DIGIT) 0;
        for (j = 0; j < DIGIT_SIZE_B; j++) {
            poly[i] |= (DIGIT) poly_bytes[i * DIGIT_SIZE_B + j] << 8 * j;
        }
    }
}

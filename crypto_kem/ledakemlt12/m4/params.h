#ifndef PARAMS_H
#define PARAMS_H

#include "fips202.h"

#include <inttypes.h>
#include <stddef.h>

// QC-LDPC parameters
#define TRNG_BYTE_LENGTH (24)
#define HASH_BYTE_LENGTH (32)
#define HASH_FUNCTION sha3_256

#define N0              (2)
#define P               (52147)  // modulus(x) = x^P-1
#define DV              (9)      // odd number
#define M               (9)
#define M0              (5)
#define M1              (4)
#define NUM_ERRORS_T    (136)

// Derived parameters, they are useful for QC-LDPC algorithms
#define HASH_BIT_LENGTH (HASH_BYTE_LENGTH << 3)
#define K               ((N0-1)*P)
#define N               (N0*P)
#define DC              (N0*DV)

#define Q_BLOCK_WEIGHTS  {{M0,M1},{M1,M0}}
static const uint8_t qBlockWeights[N0][N0] = Q_BLOCK_WEIGHTS;

/* Decoding definitions for DFR level 2^-SL with SL=128 */
#define ITERATIONS_MAX  (2)
#define B0              (43)
#define T_BAR           (4)

// Ring constants
#define NUM_BITS_GF2X_ELEMENT                   (P) // 52147
#define NUM_DIGITS_GF2X_ELEMENT                 ((P+DIGIT_SIZE_b-1)/DIGIT_SIZE_b)
#define MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT    ((P % DIGIT_SIZE_b) ? (P % DIGIT_SIZE_b)-1 : DIGIT_SIZE_b-1)
#define NUM_BITS_GF2X_MODULUS                   (P+1)
#define NUM_DIGITS_GF2X_MODULUS                 ((P+1+DIGIT_SIZE_b-1)/DIGIT_SIZE_b)
#define MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS    (P-DIGIT_SIZE_b*(NUM_DIGITS_GF2X_MODULUS-1))
#define INVALID_POS_VALUE                       (P)
#define P_BITS                                  (16) // log_2(p) = 15.6703

#define LOG_DIGIT_SIZE_b                        5
#define TAIL                                    (P % DIGIT_SIZE_b)
#define TAIL_MASK                               ((1 << TAIL) - 1)
#define RLEN                                    (POLY_DIGITS << 1)

// gf2x multiprecision constants
typedef uint32_t            DIGIT;
#define DIGIT_SIZE_B        4
#define DIGIT_SIZE_b        (DIGIT_SIZE_B << 3)
#define POSITION_T          uint32_t

// Multiplication constants for Karatsuba TC3W
#define MIN_KAR_DIGITS      10
#define MIN_TOOM_DIGITS     42
#define STACK_KAR_ONLY      4878
#define STACK_WORDS         5780

// Sizes of polynomials error polynomials 
#define POLY_DIGITS         NUM_DIGITS_GF2X_ELEMENT
#define ERROR_DIGITS        N0 * NUM_DIGITS_GF2X_ELEMENT
#define POLY_BYTES          POLY_DIGITS * DIGIT_SIZE_B
#define ERROR_BYTES         ERROR_DIGITS * DIGIT_SIZE_B

#endif

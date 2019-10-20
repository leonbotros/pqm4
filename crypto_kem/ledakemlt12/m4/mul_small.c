#include "gf2x.h"
#include "mul_small.h"

/* Small-degree multiplications
 * Implementation by Richard Brent, Pierrick Gaudry, Emmanuel Thome and Paul Zimmerman
 * Modified souce code from gf2x library for LEDAcrypt by L. Botros, August 1, 2019
 * See: gf2x.gforge.inra.fr
 */

/* mul1 in mull_small.h */

/* Memory-reduced Karatsuba variant using mul1 */
void gf2x_mul2(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT t;
    DIGIT u[2];
    gf2x_mul1(c, a[0], b[0]);
    gf2x_mul1(c + 2, a[1], b[1]);
    t = c[1] ^ c[2];
    gf2x_mul1(u, a[0] ^ a[1], b[0] ^ b[1]);
    c[1] = c[0] ^ u[0] ^ t;
    c[2] = c[3] ^ u[1] ^ t;
}

/* Weimerskirch and Paar, 2003
 * Generalizations of the Karatsuba Algorithm for Efficient Implementations
 * Uses 6 mul1 */
void gf2x_mul3(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT aa[3], bb[3];
    DIGIT p0[2], p1[2], p2[2];
    DIGIT pp0[2], pp1[2], pp2[2];
    aa[0] = a[1] ^ a[2];
    aa[1] = a[0] ^ a[2];
    aa[2] = a[0] ^ a[1];
    bb[0] = b[1] ^ b[2];
    bb[1] = b[0] ^ b[2];
    bb[2] = b[0] ^ b[1];
    gf2x_mul1(p0, a[0], b[0]);
    gf2x_mul1(p1, a[1], b[1]);
    gf2x_mul1(p2, a[2], b[2]);
    gf2x_mul1(pp0, aa[0], bb[0]);
    gf2x_mul1(pp1, aa[1], bb[1]);
    gf2x_mul1(pp2, aa[2], bb[2]);
    c[0] = p0[0];
    c[1] = p0[0]  ^ p1[0]  ^ pp2[0] ^ p0[1];
    c[2] = p0[0]  ^ p1[0]  ^ p2[0]  ^ pp1[0] ^ p0[1] ^ p1[1] ^ pp2[1];
    c[3] = pp0[0] ^ p1[0]  ^ p2[0]  ^ p0[1]  ^ p1[1] ^ p2[1] ^ pp1[1];
    c[4] = p2[0]  ^ pp0[1] ^ p1[1]  ^ p2[1];
    c[5] = p2[1];
}

/* Karatsuba with mul2 */
void gf2x_mul4(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT aa[2], bb[2], ab[4];
    DIGIT lo[4], hi[4];
    gf2x_mul2(lo, a, b);
    gf2x_mul2(hi, a + 2, b + 2);
    aa[0] = a[0] ^ a[2];
    aa[1] = a[1] ^ a[3];
    bb[0] = b[0] ^ b[2];
    bb[1] = b[1] ^ b[3];
    DIGIT c24 = lo[2] ^ hi[0];
    DIGIT c35 = lo[3] ^ hi[1];
    gf2x_mul2(ab, aa, bb);
    c[0] = lo[0];
    c[1] = lo[1];
    c[2] = ab[0] ^ lo[0] ^ c24;
    c[3] = ab[1] ^ lo[1] ^ c35;
    c[4] = ab[2] ^ hi[2] ^ c24;
    c[5] = ab[3] ^ hi[3] ^ c35;
    c[6] = hi[2];
    c[7] = hi[3];
}

/* Montgomery formulae with 13 multiplications, see
 * Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
 * IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
void gf2x_mul5(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT ta[3], tb[3], pa[8], pb[8], p[26], t[14];
    // Compute as and bs
    ta[0] = a[0]  ^ a[4];
    tb[0] = b[0]  ^ b[4];
    ta[1] = a[1]  ^ a[2];
    tb[1] = b[1]  ^ b[2];
    ta[2] = a[3]  ^ ta[0];
    tb[2] = b[3]  ^ tb[0];
    pa[0] = ta[1] ^ ta[2];
    pb[0] = tb[1] ^ tb[2];
    pa[1] = a[2]  ^ ta[2];
    pb[1] = b[2]  ^ tb[2];
    pa[2] = ta[0] ^ ta[1];
    pb[2] = tb[0] ^ tb[1];
    pa[3] = a[1]  ^ ta[2];
    pb[3] = b[1]  ^ tb[2];
    pa[4] = a[0]  ^ a[2] ^ a[3];
    pb[4] = b[0]  ^ b[2] ^ b[3];
    pa[5] = a[4]  ^ ta[1];
    pb[5] = b[4]  ^ tb[1];
    pa[6] = a[3]  ^ a[4];
    pb[6] = b[3]  ^ b[4];
    pa[7] = a[0]  ^ a[1];
    pb[7] = b[0]  ^ b[1];

    // Multiply
    gf2x_mul1(p + 0,  pa[0], pb[0]);
    gf2x_mul1(p + 2,  pa[1], pb[1]);
    gf2x_mul1(p + 4,  pa[2], pb[2]);
    gf2x_mul1(p + 6,  pa[3], pb[3]);
    gf2x_mul1(p + 8,  pa[4], pb[4]);
    gf2x_mul1(p + 10, pa[5], pb[5]);
    gf2x_mul1(p + 12, pa[6], pb[6]);
    gf2x_mul1(p + 14, pa[7], pb[7]);
    gf2x_mul1(p + 16, ta[0], tb[0]);
    gf2x_mul1(p + 18, a[4],  b[4]);
    gf2x_mul1(p + 20, a[3],  b[3]);
    gf2x_mul1(p + 22, a[1],  b[1]);
    gf2x_mul1(p + 24, a[0],  b[0]);

    t[0]  = p[14] ^ p[24];
    t[1]  = p[15] ^ p[25];
    t[2]  = p[12] ^ p[18];
    t[3]  = p[13] ^ p[19];
    t[4]  = p[2]  ^ p[16];
    t[5]  = p[3]  ^ p[17];
    t[6]  = p[0]  ^ p[6];
    t[7]  = p[1]  ^ p[7];
    t[8]  = p[4]  ^ p[16];
    t[9]  = p[5]  ^ p[17];
    t[10] = p[10] ^ t[0];
    t[11] = p[11] ^ t[1];
    t[12] = p[8]  ^ t[2];
    t[13] = p[9]  ^ t[3];

    c[0] = p[24];
    c[1] = p[22]  ^ t[0]  ^ p[25];
    c[2] = p[18]  ^ t[8]  ^ t[10]  ^ p[23] ^ t[1];
    c[3] = t[2]   ^ t[4]  ^ t[6]   ^ p[19] ^ t[9]   ^ t[11];
    c[4] = p[0]   ^ p[20] ^ p[22]  ^ t[10] ^ t[12]  ^ t[3]   ^ t[5]  ^ t[7];
    c[5] = t[0]   ^ t[6]  ^ t[8]   ^ p[1]  ^ p[21]  ^ p[23]  ^ t[11] ^ t[13];
    c[6] = p[24]  ^ t[4]  ^ t[12]  ^ t[1]  ^ t[7]   ^ t[9];
    c[7] = p[20]  ^ t[2]  ^ p[25]  ^ t[5]  ^ t[13];
    c[8] = p[18]  ^ p[21] ^ t[3];
    c[9] = p[19];
}

/* This code uses the K3 formula from Weimerskirch and Paar,
 * http://weimerskirch.org/files/Weimerskirch_Karatsuba.pdf,
 * which performs only 6 calls to gf2x_mul2.
 * Can save 1 mult using Montgomery formula */
void gf2x_mul6(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT d01[4], d1[4], d12[4], aa[2], bb[2];
    gf2x_mul2(c, a, b);                                    /* D0 */
    gf2x_mul2(d1, a + 2, b + 2);                           /* D1 */
    gf2x_mul2(c + 8, a + 4, b + 4);                        /* D2 */
    aa[0] = a[0] ^ a[2];
    aa[1] = a[1] ^ a[3];
    bb[0] = b[0] ^ b[2];
    bb[1] = b[1] ^ b[3];
    gf2x_mul2(d01, aa, bb);                                /* D01 */
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    gf2x_mul2(c + 4, aa, bb);                              /* D02 */
    aa[0] = a[2] ^ a[4];
    aa[1] = a[3] ^ a[5];
    bb[0] = b[2] ^ b[4];
    bb[1] = b[3] ^ b[5];
    gf2x_mul2(d12, aa, bb);                                /* D12 */
    c[2] ^= d1[0];
    c[3] ^= d1[1];                                         /* low(D1) + high(D0) */
    c[4] ^= c[2];
    c[5] ^= c[3];                                          /* low(D02) + low(D1) + high(D0) */
    d12[0] ^= c[2];
    d12[1] ^= c[3];                                        /* low(D12) + low(D1) + high(D0) */
    c[8] ^= d1[2];
    c[9] ^= d1[3];                                         /* low(D2) + high(D1) */
    c[6] ^= c[8];
    c[7] ^= c[9];                                          /* high(D02) + low(D2) + high(D1) */
    d01[2] ^= c[8];
    d01[3] ^= c[9];                                        /* high(D01) + low(D2) + high(D1) */
    c[2] ^= d01[0] ^ c[0];
    c[3] ^= d01[1] ^ c[1];                                 /* l(D1)+h(D0)+l(D01)+l(D0) */
    c[4] ^= c[0] ^ d01[2];
    c[5] ^= c[1] ^ d01[3];
    c[6] ^= d12[0] ^ c[10];
    c[7] ^= d12[1] ^ c[11];
    c[8] ^= d12[2] ^ c[10];
    c[9] ^= d12[3] ^ c[11];
}

/* can do save 2 mults by using seven-term formula from Montgomery */
void gf2x_mul7(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7;
    gf2x_mul3(c + 8, a + 4, b + 4);
    gf2x_mul4(c, a, b);
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3];
    gf2x_mul4(ab, aa, bb);
    ab4 = ab[4] ^ c[4];
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    c[4] ^= ab[0] ^ c[0] ^ c[8];
    c[5] ^= ab[1] ^ c[1] ^ c[9];
    c[6] ^= ab[2] ^ c[2] ^ c[10];
    c[7] ^= ab[3] ^ c[3] ^ c[11];
    c[8] ^= ab4 ^ c[12];
    c[9] ^= ab5 ^ c[13];
    c[10] ^= ab6;
    c[11] ^= ab7;
}

/* 3 calls to mul4, i.e., 27 mul1 */
void gf2x_mul8(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT aa[4], bb[4], cc[4];
    gf2x_mul4(c + 8, a + 4, b + 4);
    gf2x_mul4(c, a, b);
    cc[0] = c[4] ^ c[8];
    cc[1] = c[5] ^ c[9];
    cc[2] = c[6] ^ c[10];
    cc[3] = c[7] ^ c[11];
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3] ^ a[7];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3] ^ b[7];
    gf2x_mul4(c + 4, aa, bb);
    c[4]  ^= c[0]  ^ cc[0];
    c[5]  ^= c[1]  ^ cc[1];
    c[6]  ^= c[2]  ^ cc[2];
    c[7]  ^= c[3]  ^ cc[3];
    c[8]  ^= c[12] ^ cc[0];
    c[9]  ^= c[13] ^ cc[1];
    c[10] ^= c[14] ^ cc[2];
    c[11] ^= c[15] ^ cc[3];
}

/* 1 call to mul4 and 2 calls to mul5, i.e., 35 mul1  */
void gf2x_mul9(DIGIT *c, const DIGIT *a, const DIGIT *b) {
    DIGIT aa[5], bb[5], ab[10], ab5, ab6, ab7, ab8, ab9;
    gf2x_mul4(c + 10, a + 5, b + 5);
    gf2x_mul5(c, a, b);
    aa[0] = a[0] ^ a[5];
    aa[1] = a[1] ^ a[6];
    aa[2] = a[2] ^ a[7];
    aa[3] = a[3] ^ a[8];
    aa[4] = a[4];
    bb[0] = b[0] ^ b[5];
    bb[1] = b[1] ^ b[6];
    bb[2] = b[2] ^ b[7];
    bb[3] = b[3] ^ b[8];
    bb[4] = b[4];
    gf2x_mul5(ab, aa, bb);
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    ab8 = ab[8] ^ c[8];
    ab9 = ab[9] ^ c[9];
    c[5] ^= ab[0] ^ c[0] ^ c[10];
    c[6] ^= ab[1] ^ c[1] ^ c[11];
    c[7] ^= ab[2] ^ c[2] ^ c[12];
    c[8] ^= ab[3] ^ c[3] ^ c[13];
    c[9] ^= ab[4] ^ c[4] ^ c[14];
    c[10] ^= ab5 ^ c[15];
    c[11] ^= ab6 ^ c[16];
    c[12] ^= ab7 ^ c[17];
    c[13] ^= ab8;
    c[14] ^= ab9;
}



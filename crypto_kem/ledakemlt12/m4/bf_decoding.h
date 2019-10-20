#ifndef BF_DECODING_H
#define BF_DECODING_H

#include "gf2x.h"
#include "params.h"


int bf_decoding(DIGIT err[],
                const POSITION_T HtrPosOnes[N0][DV],
                const POSITION_T QtrPosOnes[N0][M],
                DIGIT privateSyndrome[],
                uint8_t threshold); // B2

#endif

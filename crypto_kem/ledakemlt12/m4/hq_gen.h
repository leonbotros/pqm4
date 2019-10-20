#ifndef H_Q_MATRICES_GENERATION_H
#define H_Q_MATRICES_GENERATION_H

#include "gf2x.h"
#include "params.h"
#include "rng.h"

void generateHPosOnes(POSITION_T HPosOnes[N0][DV], AES_XOF_struct *keys_expander);
void generateQPosOnes(POSITION_T QPosOnes[N0][M], AES_XOF_struct *keys_expander);
void transposeHPosOnes(POSITION_T HtrPosOnes[N0][DV], POSITION_T HPosOnes[N0][DV]);
void transposeQPosOnes(POSITION_T QtrPosOnes[N0][M], POSITION_T QPosOnes[N0][M]);

#endif

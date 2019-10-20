#include "bf_decoding.h"
#include "ring.h"

#include <string.h>

int bf_decoding(DIGIT err[],
                const POSITION_T HtrPosOnes[N0][DV],
                const POSITION_T QtrPosOnes[N0][M],
                DIGIT privateSyndrome[],
                uint8_t secondIterThreshold) {

    POSITION_T currQBlkPos[M], currQBitPos[M];
    POSITION_T syndromePosToFlip, tmp, hPos;
    uint32_t correlation, corrt_syndrome_based;
    uint8_t unsat;
    size_t currQoneIdx, endQblockIdx;
    size_t check;
    int iteration = 0;

    memset(err, 0x00, ERROR_BYTES);
    do {
        /* iteration based threshold determination*/
        corrt_syndrome_based = iteration * secondIterThreshold + (1 - iteration) * B0;

        // Computation of correlation  with a full Q matrix
        for (size_t i = 0; i < N0; i++) {
            for (POSITION_T j = 0; j < P; j++) {
                currQoneIdx =  endQblockIdx = 0;
                correlation = 0;

                for (size_t blockIdx = 0; blockIdx < N0; blockIdx++) {
                    endQblockIdx += qBlockWeights[blockIdx][i];
                    //currblockoffset = blockIdx * P;
                    for (; currQoneIdx < endQblockIdx; currQoneIdx++) {
                        tmp = QtrPosOnes[i][currQoneIdx] + j;
                        tmp = tmp >= P ? tmp - P : tmp;
                            
                        // wildy unefficient
                        unsat = 0;
                        for (size_t HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                            hPos = (HtrPosOnes[blockIdx][HtrOneIdx] + tmp) >= P ?
                                   (HtrPosOnes[blockIdx][HtrOneIdx] + tmp) - P :
                                   (HtrPosOnes[blockIdx][HtrOneIdx] + tmp);
                            if (gf2x_get_coeff(privateSyndrome, hPos)) {
                                unsat++;
                            }
                        }

                        currQBitPos[currQoneIdx] = tmp;
                        currQBlkPos[currQoneIdx] = (POSITION_T)blockIdx;
                        correlation += unsat;
                    }
                }

                /* Correlation based flipping */
                if (correlation >= corrt_syndrome_based) {
                    gf2x_toggle_coeff(err + POLY_DIGITS * i, j);
                    for (size_t v = 0; v < M; v++) {
                        for (size_t HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                            syndromePosToFlip = (HtrPosOnes[currQBlkPos[v]][HtrOneIdx] + currQBitPos[v]);
                            syndromePosToFlip = syndromePosToFlip >= P ? syndromePosToFlip - P : syndromePosToFlip;
                            gf2x_toggle_coeff(privateSyndrome, syndromePosToFlip);
                        }
                    } // end for v
                } // end if
            } // end for j
        } // end for i

        iteration = iteration + 1;
        check = 0;
        while (check < POLY_DIGITS && privateSyndrome[check++] == 0) {};

    } while (iteration < ITERATIONS_MAX && check < POLY_DIGITS);

    return (check == POLY_DIGITS);
}

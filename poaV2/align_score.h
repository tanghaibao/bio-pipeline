#ifndef ALIGN_SCORE_HEADER_INCLUDED
#define ALIGN_SCORE_HEADER_INCLUDED

#include <default.h>
#include <poa.h>
#include <seq_util.h>

/*********************************************************** align_score.c */
LPOScore_T matrix_scoring_function(int i,
				   int j,
				   LPOLetter_T seq_x[],
				   LPOLetter_T seq_y[],
				   ResidueScoreMatrix_T *m);

#endif

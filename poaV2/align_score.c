



#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"



typedef struct {
  unsigned char x;
  unsigned char y:7;
  unsigned char is_aligned:1;
} LPOMove_T;



/* YOU CAN PUT ANY SCORING METHOD YOU WANT INSIDE THIS
   FUNCTION. JUST REPLACE THE CONTENTS OF THE FUNCTION WITH
   YOUR SCORING METHOD */
LPOScore_T matrix_scoring_function(int i,
				   int j,
				   LPOLetter_T seq_x[],
				   LPOLetter_T seq_y[],
				   ResidueScoreMatrix_T *m)
{
  return m->score[seq_x[i].letter][seq_y[j].letter]; /*TRIVIAL SCORING FUNC:
						       JUST USE MATRIX VALUE*/
}

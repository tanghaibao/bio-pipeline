


#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"





typedef struct {
  unsigned char x;
  unsigned char y:7;
  unsigned char is_aligned:1;
} LPOMove_T;


typedef unsigned char LPOGapLength_T;





void trace_back_lpo_po_alignment(int len_x,LPOLetter_T seq_x[],
				 int len_y,LPOLetter_T seq_y[],
				 LPOMove_T **move,
				 LPOLetterRef_T best_x,
				 LPOLetterRef_T best_y,
				 LPOLetterRef_T **x_to_y,
				 LPOLetterRef_T **y_to_x)
{
  int i;
  LPOLetterRef_T new_x,*x_al=NULL,*y_al=NULL;
  LPOLetterLink_T *left;

  
  CALLOC(x_al,len_x,LPOLetterRef_T);
  CALLOC(y_al,len_y,LPOLetterRef_T);
  LOOP (i,len_x) x_al[i]= INVALID_LETTER_POSITION;
  LOOP (i,len_y) y_al[i]= INVALID_LETTER_POSITION;

  while (best_x>=0 && best_y>=0) {
    if (move[best_x][best_y].is_aligned) {/* ALIGNED! MAP best_x <--> best_y */
      x_al[best_x]=best_y;
      y_al[best_y]=best_x;
    }

    if (0==move[best_x][best_y].x /* HIT START OF THE ALIGNMENT, SO QUIT */
	&& 0==move[best_x][best_y].y)
      break;

    if ((i=move[best_x][best_y].x)>0) { /* TRACE BACK ON X */
      for (left= &seq_x[best_x].left;--i >0;left=left->more); /* USE iTH MOVE*/
      new_x = left->ipos;
    }
    else new_x=best_x;
    if ((i=move[best_x][best_y].y)>0) { /* TRACE BACK ON Y */
      for (left= &seq_y[best_y].left;--i >0;left=left->more); /* USE iTH MOVE*/
      best_y = left->ipos;
    }
    best_x=new_x;
  }


  if (x_to_y) /* HAND BACK ALIGNMENT RECIPROCAL MAPPINGS */
    *x_to_y = x_al;
  else
    free(x_al);
  if (y_to_x)
    *y_to_x = y_al;
  else
    free(y_al);
  return;
}



/** performs partial order alignment:
  seq_x[] may be a partial order;
  seq_y[] may be a partial order;
  returns the alignment in x_to_y[] and y_to_x, and also 
  returns the alignment score as the return value */
LPOScore_T align_lpo_po (LPOSequence_T *lposeq_x,
			 LPOSequence_T *lposeq_y,
			 ResidueScoreMatrix_T *m,
			 LPOLetterRef_T **x_to_y,
			 LPOLetterRef_T **y_to_x,
			 LPOScore_T (*scoring_function)
			 (int,int,LPOLetter_T [],LPOLetter_T [],ResidueScoreMatrix_T *),
			 int use_global_alignment)
{
  int len_x = lposeq_x->length;
  int len_y = lposeq_y->length;
  LPOLetter_T *seq_x = lposeq_x->letter;
  LPOLetter_T *seq_y = lposeq_y->letter;
  
  int i,j,j_left,best_x,best_y,nfree_list=0,ilast_right=0;
  LPOScore_T **score=NULL,*my_score,*my_matrix,*score_base=NULL;
  LPOMove_T **move=NULL,*move_base=NULL,*my_move;
  LPOGapLength_T **gap_length=NULL,*my_gap_length,*gap_length_base=NULL;
  LPOLetterLink_T *left,*my_left,*y_left;
  int new_gap_len,insert_x,previous_x,previous_y,
    i_x,new_score,new_x,new_y,current_gap_length,i_y,insert_y;
  LPOScore_T match_score,previous_score,
    insert_x_try,insert_x_score,insert_y_score,best_score= -999999;

  int max_gap_length;
  LPOScore_T *gap_penalty_x, *gap_penalty_y;  
  
  max_gap_length = m->max_gap_length;
  gap_penalty_x = m->gap_penalty_x;
  gap_penalty_y = m->gap_penalty_y;
  
  CALLOC(score,len_x,LPOScore_T *); /* ALLOCATE MATRIX STORAGE: ROW POINTERS */
  CALLOC(move,len_x,LPOMove_T *);
  CALLOC(gap_length,len_x,LPOGapLength_T *);
  CALLOC(score_base,len_x*(len_y+1),LPOScore_T); /*ALLOCATE MATRIX RECTANGLE */
  CALLOC(move_base,len_x*(len_y+1),LPOMove_T); /*ALLOCATE MATRIX RECTANGLE */
  CALLOC(gap_length_base,len_x*(len_y+1),LPOGapLength_T); /*ALLOCATE MATRIX RECTANGLE */


  LOOPF (i,len_x) {/* BUILD UP DP MATRIX, ROW BY ROW */
    score[i]=score_base+i*(len_y+1)+1; /* LEAVE SPACE FOR [-1] ENTRY */
    move[i]=move_base+i*(len_y+1)+1; /* LEAVE SPACE FOR [-1] ENTRY */
    gap_length[i]=gap_length_base+i*(len_y+1)+1; /*LEAVE SPACE FOR [-1] ENTRY*/
    my_move=move[i]; /* USED TO SPEED UP MATRIX ACCESS INSIDE INNER LOOP*/
    my_score=score[i];
    my_gap_length=gap_length[i];
    score[i][-1]= -999999; /* UNACCEPTABLE SCORE ENSURES -1 NEVER CHOSEN*/
    /* my_matrix=m->score[seq_x[i].letter]; NOT USED */
    if (seq_x[i].left.ipos>=0) /* AT LEAST ONE VALID POSITION TO THE LEFT */
      my_left= &seq_x[i].left;
    else /* THERE IS NO POSITION TO THE LEFT */
      my_left=NULL;
    LOOPF (j,len_y) {
      j_left=SEQ_Y_LEFT(j); /* POSITION TO THE LEFT OF j */
      previous_score=previous_x=previous_y=0;
      i_x=1;
      insert_x_score= -999999;
      for (left=my_left;left;left=left->more) {
	if (move[left->ipos][j].x>0) /* COULD BE [X,0] GAP CONTINUATION */
	  current_gap_length=gap_length[left->ipos][j];
	else/*NOT AN EXTENSION OF A [X,0] GAP, SO TREAT AS START OF NEW GAP*/
	  current_gap_length=0;
	insert_x_try=score[left->ipos][j] /* FIND BEST [X,0] MOVE */
	  + left->score /* INCLUDE WEIGHTING FROM THIS EDGE */
	  - gap_penalty_x[current_gap_length];
	if (insert_x_try>insert_x_score) { /* IF BEST insert_x MOVE, SAVE*/
	  insert_x=i_x;
	  insert_x_score=insert_x_try;
	  new_gap_len=current_gap_length+1;
	  if (new_gap_len>max_gap_length) /* PREVENT OVERFLOW */
	    new_gap_len=max_gap_length;
	}

	 /* FIND BEST [X,Y] MOVE */
	if (seq_y[j].left.ipos>=0){/*AT LEAST ONE VALID POSITION TO THE LEFT*/
	  i_y=1;
	  for (y_left= &seq_y[j].left;y_left;y_left=y_left->more) {
	    if (score[left->ipos][y_left->ipos] + left->score + y_left->score
		>previous_score) {
	      previous_score=score[left->ipos][y_left->ipos] 
		+ left->score + y_left->score;
	      previous_x=i_x;
	      previous_y=i_y;
	    }
	    i_y++;
	  }
	}
	i_x++; /* ADVANCE X MOVE INDEX */
      } /* DONE SCANNING PREDECESSORS ON X */

      match_score = previous_score /* TAKE BEST PREDECESSOR */
	+ scoring_function(i,j,seq_x,seq_y,m);
#ifdef USE_LOCAL_NEUTRALITY_CORRECTION /* NO LONGER USED */
      if (seq_x[i].score<seq_y[j].score) /* USE STRONGEST NEUTRALITY ADJ'T*/
	match_score += seq_x[i].score;
      else
	match_score += seq_y[j].score;
#endif
      if (match_score>insert_x_score) { /* PREFER [X,Y] MOVE */
	new_score=match_score;
	new_x=previous_x;
	new_y=previous_y;
	new_gap_len=0;
      }
      else { /* PREFER [X,0] MOVE */
	new_score=insert_x_score;
	new_x=insert_x;
	new_y=0;
      }

       /* [0,Y] MOVE */
      insert_y_score= -999999;
      if (seq_y[j].left.ipos>=0){/*AT LEAST ONE VALID POSITION TO THE LEFT*/
	i_y=1;
	for (y_left= &seq_y[j].left;y_left;y_left=y_left->more) {
	  if (my_move[y_left->ipos].y>0) /* COULD BE [0,1] GAP CONTINUATION */
	    current_gap_length=my_gap_length[y_left->ipos];
	  else/*NOT AN EXTENSION OF A [0,1] GAP, SO TREAT AS START OF NEW GAP*/
	    current_gap_length=0;
	  if (insert_y_score <
	      my_score[y_left->ipos]-gap_penalty_y[current_gap_length]) {
	    insert_y_score=my_score[y_left->ipos]-gap_penalty_y[current_gap_length];
	    insert_y=i_y;
	  }
	  i_y++;
	}
      }

      if (insert_y_score<new_score) { /* [new_x,new_y] MOVE IS BEST */
	my_score[j]=new_score;
	my_move[j].x=new_x;
	my_move[j].y=new_y;
	my_gap_length[j]=new_gap_len;
	if (new_gap_len==0)
	  my_move[j].is_aligned=1;
      }
      else { /* [0,Y] MOVE IS BEST */
	my_score[j]=insert_y_score;
	my_move[j].x=0;
	my_move[j].y=insert_y;
	my_gap_length[j]=current_gap_length+1;
	if (my_gap_length[j]>max_gap_length) /* PREVENT OVERFLOW */
	  my_gap_length[j]=max_gap_length;
      }
      if (my_score[j]>best_score) { /* RECORD BEST MOVE */
	best_score=my_score[j];
	best_x=i;
	best_y=j;
      }
    }
  } /* DYNAMIC PROGRAMING MATRIX COMPLETE, NOW TRACE BACK FROM best_x,best_y*/

  IF_GUARD(best_x>=len_x || best_y>=len_y,1.1,(ERRTXT,"Bounds exceeded!\nbest_x,best_y:%d,%d\tlen:%d,%d\n",best_x,best_y,len_x,len_y),CRASH);
  trace_back_lpo_po_alignment(len_x,seq_x,len_y,seq_y,move,best_x,best_y,
			      x_to_y,y_to_x);
  

  FREE(score_base); /* FREE MEMORY */
  FREE(move_base);
  FREE(gap_length_base);
  FREE(score);
  FREE(move);
  FREE(gap_length);
  return best_score;
}






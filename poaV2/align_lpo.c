

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"



typedef struct {
  unsigned char x:6;
  unsigned char y:1;
  unsigned char is_aligned:1;
} LPOMove_T;


typedef unsigned char LPOGapLength_T;





void trace_back_lpo_alignment(int len_x,LPOLetter_T seq_x[],
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
      for (left= &seq_x[best_x].left;i-- >0;left=left->more) /* USE iTH MOVE*/
	new_x = left->ipos;
    }
    else /* NO MOVE ON X */
      new_x=best_x;

    if (move[best_x][best_y].y>0) /* ASSUMING seq_y A LINEAR SEQUENCE*/
      best_y=SEQ_Y_LEFT(best_y); /* TRACE BACK ON Y */
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





/* FOR CONTROLLING THE FREEING OF score[] ARRAYS */
typedef struct {
  LPOLetterRef_T last_right; /* LAST POSITION WHICH REFERENCES THIS POSITION*/
  LPOLetterRef_T my_pos; /* INDEX OF THIS POSITION */
} LastRightList_T;


/* SORT IN ASCENDING ORDER BY last_right */
int last_right_qsort_cmp(const void *void_a,const void *void_b)
{
  const LastRightList_T *a=(const LastRightList_T *)void_a,
  *b=(const LastRightList_T *)void_b;
  if (a->last_right < b->last_right)
    return -1;
  if (a->last_right == b->last_right)
    return 0;
  else
    return 1;
}




LastRightList_T *last_right_list(int len,LPOLetter_T seq[])
{
  int i;
  LastRightList_T *list=NULL;
  LPOLetterLink_T *right;

  CALLOC(list,len,LastRightList_T);
  LOOP (i,len) {
    list[i].last_right=list[i].my_pos=i; /* DEFAULT: NOTHING RIGHT OF THIS*/
    for (right= &seq[i].right;right && right->ipos>=0;right=right->more)
      if (right->ipos>list[i].last_right)
	list[i].last_right=right->ipos;
  }
  qsort(list,len,sizeof(LastRightList_T),last_right_qsort_cmp);
  return list;
}





typedef struct {
  LPOScore_T *score;
  LPOGapLength_T *gap_length;
} LPORowFreeList_T;


void free_row_list(int nfree_list,LPORowFreeList_T free_list[])
{
/*  fprintf(stderr,"Maximum #rows allocated was %d\n\n",nfree_list);*/
  while (nfree_list-- >0) { /* DUMP EVERYTHING STORED ON THE FREE LIST */
    free(free_list[nfree_list].score);
    free(free_list[nfree_list].gap_length);
  }
}





/** performs partial order alignment:
  seq_x[] may be a partial order;
  seq_y[] is assumed to be a linear order (regular sequence);
  returns the alignment in x_to_y[] and y_to_x, and also 
  returns the alignment score as the return value */
LPOScore_T align_lpo (LPOSequence_T *lposeq_x,
		      LPOSequence_T *lposeq_y,
		      ResidueScoreMatrix_T *m,
		      LPOLetterRef_T **x_to_y,
		      LPOLetterRef_T **y_to_x,
		      int use_global_alignment)
{
  int len_x = lposeq_x->length;
  int len_y = lposeq_y->length;
  LPOLetter_T *seq_x = lposeq_x->letter;
  LPOLetter_T *seq_y = lposeq_y->letter;

  int i,j,j_left,best_x,best_y,nfree_list=0,ilast_right=0;
  LPOScore_T **score=NULL,*my_score,*my_matrix;
  LPOMove_T **move=NULL,*move_base=NULL,*my_move;
  LPOGapLength_T **gap_length=NULL,*my_gap_length;
  LPOLetterLink_T *left,*my_left;
  int new_gap_len,insert_x,previous_x,previous_y,
    i_x,new_score,new_x,new_y,current_gap_length;
  LPOScore_T match_score,previous_score,
    insert_x_try,insert_x_score,insert_y_score,best_score= -999999;
  LPORowFreeList_T *free_list=NULL;
  LastRightList_T *last_right=NULL;

  int max_gap_length;
  LPOScore_T *gap_penalty_x, *gap_penalty_y;  
  
  max_gap_length = m->max_gap_length;
  gap_penalty_x = m->gap_penalty_x;
  gap_penalty_y = m->gap_penalty_y;
  
  CALLOC(score,len_x,LPOScore_T *); /* ALLOCATE MATRIX STORAGE: ROW POINTERS */
  CALLOC(move,len_x,LPOMove_T *);
  CALLOC(gap_length,len_x,LPOGapLength_T *);
  CALLOC(free_list,len_x,LPORowFreeList_T);
  CALLOC(move_base,len_x*(len_y+1),LPOMove_T); /*ALLOCATE MATRIX RECTANGLE */
  last_right=last_right_list(len_x,seq_x); /* GET SORTED LIST TO CONTROL FREE*/

  LOOPF (i,len_x) {/* BUILD UP DP MATRIX, ROW BY ROW */
    if (nfree_list>0) { /* TAKE THE NEW ROW FROM THE FREE LIST */
      nfree_list--; /* MOVE BACK TO LAST FREE LIST ENTRY */
      score[i]=free_list[nfree_list].score+1; /* LEAVE SPACE FOR [-1] ENTRY */
      gap_length[i]=free_list[nfree_list].gap_length+1;
    }
    else { /* NEED TO ALLOCATE A NEW ROW */
      CALLOC(score[i],len_y+1,LPOScore_T);
      score[i]++; /* LEAVE SPACE FOR [-1] ENTRY */
      CALLOC(gap_length[i],len_y+1,LPOGapLength_T);
      gap_length[i]++; /* LEAVE SPACE FOR [-1] ENTRY */
    }
    move[i]=move_base+i*(len_y+1)+1; /* LEAVE SPACE FOR [-1] ENTRY */
    my_move=move[i]; /* USED TO SPEED UP MATRIX ACCESS INSIDE INNER LOOP*/
    my_score=score[i];
    my_gap_length=gap_length[i];
    score[i][-1]= -999999; /* UNACCEPTABLE SCORE ENSURES -1 NEVER CHOSEN*/
    my_matrix=m->score[seq_x[i].letter];
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
	  if (new_gap_len>max_gap_length)  /* PREVENT OVERFLOW */
	    new_gap_len=max_gap_length;
	}

	 /* FIND BEST [X,1] MOVE */
	if (score[left->ipos][j_left]+left->score>previous_score) {
	  previous_score=score[left->ipos][j_left]+left->score;
	  previous_x=i_x;
	  previous_y=1; /* ASSUMING seq_y JUST LINEAR SEQUENCE */
	}
	i_x++; /* ADVANCE X MOVE INDEX */
      } /* DONE SCANNING PREDECESSORS ON X */

      match_score = previous_score /* TAKE BEST PREDECESSOR */
	+ my_matrix[seq_y[j].letter];

      if (match_score>insert_x_score) { /* PREFER [X,1] MOVE */
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

       /* [0,1] MOVE */
      if (my_move[j_left].y==1) /* COULD BE [0,1] GAP CONTINUATION */
	current_gap_length=my_gap_length[j_left];
      else /* NOT AN EXTENSION OF A [0,1] GAP, SO TREAT AS START OF NEW GAP*/
	current_gap_length=0;
      insert_y_score=my_score[j_left]-gap_penalty_y[current_gap_length];

      if (insert_y_score<new_score) { /* [new_x,new_y] MOVE IS BEST */
	my_score[j]=new_score;
	my_move[j].x=new_x;
	my_move[j].y=new_y;
	my_gap_length[j]=new_gap_len;
	if (new_gap_len==0)
	  my_move[j].is_aligned=1;
      }
      else { /* [0,1] MOVE IS BEST */
	my_score[j]=insert_y_score;
	my_move[j].x=0;
	my_move[j].y=1;
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
    while (ilast_right<len_x && last_right[ilast_right].last_right<=i) {
      free_list[nfree_list].gap_length= /* PUSH ROW ONTO FREE LIST FOR REUSE*/
	gap_length[last_right[ilast_right].my_pos]-1;/*READJUST TO BASE ADDR!*/
      free_list[nfree_list].score=score[last_right[ilast_right].my_pos]-1;
      score[last_right[ilast_right].my_pos]=NULL; /*PREVENT DANGLING POINTER!*/
      gap_length[last_right[ilast_right].my_pos]=NULL;
      nfree_list++; /* INCREMENT FREE LIST SIZE */
      ilast_right++; /* MOVE TO THE NEXT POTENTIAL ROW FOR FREEING */
    }
  } /* DYNAMIC PROGRAMING MATRIX COMPLETE, NOW TRACE BACK FROM best_x,best_y*/

  IF_GUARD(best_x>=len_x || best_y>=len_y,1.1,(ERRTXT,"Bounds exceeded!\nbest_x,best_y:%d,%d\tlen:%d,%d\n",best_x,best_y,len_x,len_y),CRASH);
  trace_back_lpo_alignment(len_x,seq_x,len_y,seq_y,move,best_x,best_y,
			   x_to_y,y_to_x);


  free_row_list(nfree_list,free_list);
  FREE(move_base); /* DUMP ALLOCATED MATRIX */
  FREE(score);
  FREE(move);
  FREE(gap_length);
  FREE(free_list);
  FREE(last_right);
  return best_score;
}






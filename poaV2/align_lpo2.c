
#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/** set nonzero for old scoring (gap-opening penalty for X-Y transition) */
#define DOUBLE_GAP_SCORING (0)


typedef struct {
  unsigned char x:7;
  unsigned char y:1;
}
DPMove_T;

typedef struct {
  LPOScore_T score;
  short gap_x, gap_y;
}
DPScore_T;

/** is node 'i' the first node in any sequence in lposeq? */
static int is_initial_node (int i, LPOSequence_T *lposeq)
{
  LPOLetterSource_T *src = &((lposeq->letter[i]).source);
  while (src != NULL && src->iseq >= 0) {
    if (src->ipos == 0) {
      return 1;
    }
    src = src->more;
  }
  return 0;
}

/** is node 'i' the last node in any sequence in lposeq? */
static int is_final_node (int i, LPOSequence_T *lposeq)
{
  LPOLetterSource_T *src = &((lposeq->letter[i]).source);
  while (src != NULL && src->iseq >= 0) {
    if (src->ipos == (lposeq->source_seq[src->iseq]).length - 1) {
      return 1;
    }
    src = src->more;
  }
  return 0;
}


static void get_seq_left_and_final (LPOSequence_T *lposeq_x,
				    LPOLetterLink_T ***x_left_ptr,
				    int **is_final_node_ptr)
{
  int i, len_x = lposeq_x->length;
  LPOLetter_T *seq_x = lposeq_x->letter;
  int *is_final_node_x = NULL;
  LPOLetterLink_T **x_left = NULL;
  
  CALLOC (is_final_node_x, len_x, int);
  
  for (i=0; i<len_x; i++) {
    is_final_node_x[i] = is_final_node (i, lposeq_x);
  }
  
  CALLOC (x_left, len_x, LPOLetterLink_T *);
  
  for (i=0; i<len_x; i++) {
    if (is_initial_node (i, lposeq_x) && seq_x[i].left.ipos != -1) {
      CALLOC (x_left[i], 1, LPOLetterLink_T);
      x_left[i]->ipos = -1;
      x_left[i]->score = 0;
      x_left[i]->more = &seq_x[i].left;
    }
    else {
      x_left[i] = &seq_x[i].left;
    }
  }

  (*is_final_node_ptr) = is_final_node_x;
  (*x_left_ptr) = x_left;
}

  
static void trace_back_lpo_alignment (int len_x, int len_y,
				      DPMove_T **move,
				      LPOLetterLink_T **x_left,
				      LPOLetterRef_T best_x, LPOLetterRef_T best_y,
				      LPOLetterRef_T **x_to_y,
				      LPOLetterRef_T **y_to_x)
{
  int i, xmove, ymove;
  LPOLetterRef_T *x_al = NULL, *y_al = NULL;
  LPOLetterLink_T *left;
  
  CALLOC (x_al, len_x, LPOLetterRef_T);
  CALLOC (y_al, len_y, LPOLetterRef_T);
  LOOP (i,len_x) x_al[i] = INVALID_LETTER_POSITION;
  LOOP (i,len_y) y_al[i] = INVALID_LETTER_POSITION;
  
  while (best_x >= 0 && best_y >= 0) {

    xmove = move[best_y][best_x].x;
    ymove = move[best_y][best_x].y;
    
    if (xmove>0 && ymove>0) { /* ALIGNED! MAP best_x <--> best_y */
      x_al[best_x]=best_y;
      y_al[best_y]=best_x;
    }

    if (xmove == 0 && ymove == 0) { /* FIRST ALIGNED PAIR */
      x_al[best_x]=best_y;
      y_al[best_y]=best_x;
      break;  /* FOUND START OF ALIGNED REGION, SO WE'RE DONE */
    }

    if (xmove>0) { /* TRACE BACK ON X */
      left = x_left[best_x];
      while ((--xmove)>0) {
	left = left->more;
      }
      best_x = left->ipos;
    }
    
    if (ymove>0) { /* TRACE BACK ON Y */
      best_y--;
    }
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
    lposeq_x may be a partial order;
    lposeq_y is assumed to be a linear order (regular sequence);
    returns the alignment in x_to_y[] and y_to_x[], and also 
    returns the alignment score as the return value.
*/

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
  
  int i, j, xcount, prev_gap, next_gap;
  int best_x = -2, best_y = -2;
  LPOScore_T best_score = -999999;
  int *is_final_node_x;
  LPOLetterLink_T **x_left = NULL, *xl;
  DPMove_T **move = NULL, *my_move;
  
  DPScore_T *curr_score = NULL, *prev_score = NULL, *init_col_score = NULL, *my_score, *swap;
  
  int max_gap_length;
  LPOScore_T *gap_penalty_x, *gap_penalty_y;
  int *next_gap_array, *next_perp_gap_array;
  
  LPOScore_T try_score, insert_x_score, insert_y_score, match_score;
  int insert_x_x, insert_x_gap;
  int insert_y_y, insert_y_gap;
  int match_x, match_y;
  
  long n_edges = 0;
  
  max_gap_length = m->max_gap_length;
  gap_penalty_x = m->gap_penalty_x;
  gap_penalty_y = m->gap_penalty_y;
  CALLOC (next_gap_array, max_gap_length + 2, int);
  CALLOC (next_perp_gap_array, max_gap_length + 2, int);

  /* GAP LENGTH EXTENSION RULE: */
  /* 0->1, 1->2, 2->3, ..., M-1->M, M->M, M+1->M+1. */
  for (i=0; i<max_gap_length+1; i++) {
    next_gap_array[i] = (i<max_gap_length) ? i+1 : i;
    next_perp_gap_array[i] = (DOUBLE_GAP_SCORING ? 0 : next_gap_array[i]);
  }
  next_gap_array[max_gap_length+1] = max_gap_length+1;
  next_perp_gap_array[max_gap_length+1] = max_gap_length+1;
  
  get_seq_left_and_final (lposeq_x, &x_left, &is_final_node_x);
  
  CALLOC (move, len_y, DPMove_T *);
  for (i=0; i<len_y; i++) {
    CALLOC (move[i], len_x, DPMove_T);
  }
  
  CALLOC (init_col_score, len_y+1, DPScore_T);
  init_col_score = &(init_col_score[1]);
  
  CALLOC (curr_score, len_x+1, DPScore_T);
  curr_score = &(curr_score[1]);
  CALLOC (prev_score, len_x+1, DPScore_T);
  prev_score = &(prev_score[1]);
  
  /* FILL INITIAL ROW. */
  /* GLOBAL ALIGNMENT: no free gaps */
  /* LOCAL ALIGNMENT: free initial gap */
  
  curr_score[-1].score = 0;
  curr_score[-1].gap_x = curr_score[-1].gap_y = 
    (use_global_alignment) ? 0 : TRUNCATE_GAP_LENGTH+1;
  
  for (i=0; i<len_x; i++) {
    for (xcount = 1, xl = x_left[i]; xl != NULL; xcount++, xl = xl->more) {
      prev_gap = curr_score[xl->ipos].gap_x;
      try_score = curr_score[xl->ipos].score + xl->score - gap_penalty_x[prev_gap];
      if (xcount == 1 || try_score > curr_score[i].score) {
	curr_score[i].score = try_score;
	curr_score[i].gap_x = next_gap_array[prev_gap];
	curr_score[i].gap_y = next_perp_gap_array[prev_gap];
      }
    }
  }
  
  /* FILL INITIAL COLUMN. */
  
  init_col_score[-1] = curr_score[-1];
  for (i=0; i<len_y; i++) {
    prev_gap = init_col_score[i-1].gap_y;
    init_col_score[i].score = init_col_score[i-1].score - gap_penalty_y[prev_gap];
    init_col_score[i].gap_x = next_perp_gap_array[prev_gap];
    init_col_score[i].gap_y = next_gap_array[prev_gap];
  }
  
  
  /** MAIN DYNAMIC PROGRAMMING LOOP **/


  /* OUTER LOOP (i-th position in linear seq y): */
  for (i=0; i<len_y; i++) {
    
    swap = prev_score; prev_score = curr_score; curr_score = swap;
    
    curr_score[-1] = init_col_score[i];
    
    /* INNER LOOP (j-th position in LPO x): */
    for (j=0; j<len_x; j++) {
      
      /* ONLY POSSIBLE Y-INSERTION: trace back to (i-1, j). */
      prev_gap = prev_score[j].gap_y;
      insert_y_score = prev_score[j].score - gap_penalty_y[prev_gap];
      insert_y_gap = prev_gap;
      
      /* ONLY ONE Y-PREDECESSOR */
      match_y = insert_y_y = 1;

      /* LOOP OVER x-predecessors: */
      for (xcount = 1, xl = x_left[j]; xl != NULL; xcount++, xl = xl->more) {
	
	/* IMPROVE XY-MATCH?: trace back to (i-1, j'=xl->ipos) */
	try_score = prev_score[xl->ipos].score + xl->score;
	if (xcount == 1 || try_score > match_score) {
	  match_score = try_score;
	  match_x = xcount;
	}
	
	/* IMPROVE X-INSERTION?: trace back to (i, j'=xl->ipos) */
	prev_gap = curr_score[xl->ipos].gap_x;
	try_score = curr_score[xl->ipos].score + xl->score - gap_penalty_x[prev_gap];
	if (xcount == 1 || try_score > insert_x_score) {
	  insert_x_score = try_score;
	  insert_x_x = xcount;
	  insert_x_gap = prev_gap;
	}
      }

      if (0 == use_global_alignment && match_score <= 0) {
	match_score = 0;
	match_x = match_y = 0;  /* FIRST ALIGNED PAIR */
      }
      
      n_edges += (xcount-1);
      match_score += m->score[(int)seq_x[j].letter][(int)seq_y[i].letter];
      
      my_score = &curr_score[j];
      my_move = &move[i][j];
      
      if (match_score > insert_y_score && match_score > insert_x_score) {
	/* XY-MATCH */
	my_score->score = match_score;
	my_score->gap_x = 0;
	my_score->gap_y = 0;
	my_move->x = match_x;
	my_move->y = match_y;
      }
      else if (insert_x_score > insert_y_score) {
	/* X-INSERTION */
	my_score->score = insert_x_score;
	my_score->gap_x = next_gap_array[insert_x_gap];
	my_score->gap_y = next_perp_gap_array[insert_x_gap];
	my_move->x = insert_x_x;
	my_move->y = 0;
      }
      else {
	/* Y-INSERTION */
	my_score->score = insert_y_score;
	my_score->gap_x = next_perp_gap_array[insert_y_gap];
	my_score->gap_y = next_gap_array[insert_y_gap];
	my_move->x = 0;
	my_move->y = insert_y_y;
      }

      /* RECORD BEST START FOR TRACEBACK */
      /* KEEPING ONLY FINAL-FINAL BESTS FOR GLOBAL ALIGNMENT */    
      if (my_score->score >= best_score && (0 == use_global_alignment || (is_final_node_x[j] && i==len_y-1))) {
	if (my_score->score > best_score || (j == best_x && i < best_y) || j < best_x) {
	  best_score = my_score->score;
	  best_x = j;
	  best_y = i;
	}
      }
    }
  }
  
  IF_GUARD(best_x>=len_x || best_y>=len_y,1.1,(ERRTXT,"Bounds exceeded!\nbest_x,best_y:%d,%d\tlen:%d,%d\n",best_x,best_y,len_x,len_y),CRASH);
  
  /*
    fprintf (stderr, "aligned (%d nodes, %ld edges) to (%d nodes)\n", len_x, n_edges/len_y, len_y);
    fprintf (stderr, "best score %d @ (%d %d)\n", best_score, best_x, best_y);
  */
  
  /* DYNAMIC PROGRAMING MATRIX COMPLETE, NOW TRACE BACK FROM best_x, best_y */
  trace_back_lpo_alignment (len_x, len_y, move, x_left,
			    best_x, best_y,
			    x_to_y, y_to_x);
  
  
  FREE (next_gap_array);
  FREE (next_perp_gap_array);
  
  prev_score = &(prev_score[-1]);
  FREE (prev_score);
  curr_score = &(curr_score[-1]);
  FREE (curr_score);
  
  init_col_score = &(init_col_score[-1]);
  FREE (init_col_score);
  
  FREE (is_final_node_x);
  
  for (i=0; i<len_x; i++) {
    if (x_left[i] != NULL && x_left[i] != &seq_x[i].left) {
      FREE (x_left[i]);
    }
  }
  FREE (x_left);
  
  for (i=0; i<len_y; i++) {
    FREE (move[i]);
  }
  FREE (move);
  
  return best_score;
}

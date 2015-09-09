
#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/** set nonzero for old scoring (gap-opening penalty for X-Y transition) */
#define DOUBLE_GAP_SCORING (0)


typedef struct {
  unsigned char x;
  unsigned char y;
}
DPMove_T;

typedef struct {
  LPOScore_T score;
  short gap_x, gap_y;
}
DPScore_T;



#define LPO_INITIAL_NODE 1
#define LPO_FINAL_NODE 2

static void get_lpo_stats (LPOSequence_T *lposeq,
			   int *n_nodes_ptr, int *n_edges_ptr, int **node_type_ptr,
			   int **refs_from_right_ptr, int *max_rows_alloced_ptr,
			   LPOLetterLink_T ***left_links_ptr)
{
  int i, j, rows_alloced = 0, max_rows_alloced = 0, n_edges = 0, len = lposeq->length;
  LPOLetter_T *seq = lposeq->letter;
  int *node_type;
  int *refs_from_right, *tmp;
  LPOLetterSource_T *src;
  LPOLetterLink_T *lnk, **left_links;
  
  CALLOC (node_type, len, int);
  CALLOC (refs_from_right, len, int);
  CALLOC (tmp, len, int);
  CALLOC (left_links, len, LPOLetterLink_T *);
    
  for (i=0; i<len; i++) {

    /* NODES CONTAINING THE FIRST RESIDUE IN ANY SEQ ARE 'INITIAL'; */
    /* DITTO, LAST RESIDUE IN ANY SEQ, 'FINAL'. */
    for (src = &(seq[i].source); src != NULL && src->iseq >= 0; src = src->more) {
      if (src->ipos == 0) {
	node_type[i] = (node_type[i] | LPO_INITIAL_NODE);
      }
      if (src->ipos == (lposeq->source_seq[src->iseq]).length - 1) {
	node_type[i] = (node_type[i] | LPO_FINAL_NODE);
      }
    }

    /* COUNTING THE LEFT-LINKS BACK TO EACH NODE ALLOWS FOR EFFICIENT */
    /* MEMORY MANAGEMENT OF 'SCORE' ROWS (in align_lpo_po). */
    for (lnk = &(seq[i].left); lnk != NULL && lnk->ipos >= 0; lnk = lnk->more) {
      refs_from_right[lnk->ipos]++;
      n_edges++;
    }
  }

  /* ALL 'INITIAL' NODES (1st in some seq) MUST BE LEFT-LINKED TO -1. */
  /* THIS ALLOWS FREE ALIGNMENT TO ANY 'BRANCH' IN GLOBAL ALIGNMENT. */
  for (i=0; i<len; i++) {
    if ((node_type[i] & LPO_INITIAL_NODE) && seq[i].left.ipos != -1) {
      CALLOC (left_links[i], 1, LPOLetterLink_T);
      left_links[i]->ipos = -1;
      left_links[i]->score = 0;
      left_links[i]->more = &seq[i].left;
    }
    else {
      left_links[i] = &seq[i].left;
    }
  }
  
  for (i=0; i<len; i++) {
    tmp[i] = refs_from_right[i];
  }
  
  for (i=0; i<len; i++) {
    rows_alloced++;
    if (rows_alloced > max_rows_alloced) {
      max_rows_alloced = rows_alloced;
    }
    for (lnk = &(seq[i].left); lnk != NULL && lnk->ipos >= 0; lnk = lnk->more) {
      if ((--tmp[lnk->ipos]) == 0) {
	rows_alloced--;
      }
    }
  }

  FREE (tmp);

  (*n_nodes_ptr) = len;
  (*n_edges_ptr) = n_edges;
  (*node_type_ptr) = node_type;
  (*refs_from_right_ptr) = refs_from_right;
  (*max_rows_alloced_ptr) = max_rows_alloced;
  (*left_links_ptr) = left_links;
}


static void trace_back_lpo_alignment (int len_x, int len_y,
				      DPMove_T **move,
				      LPOLetterLink_T **x_left,
				      LPOLetterLink_T **y_left,
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
      left = y_left[best_y];
      while ((--ymove)>0) {
	left = left->more;
      }
      best_y = left->ipos;
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


/** (align_lpo_po:)
    performs partial order alignment:
    lposeq_x and lposeq_y are partial orders;
    returns the alignment in x_to_y[] and y_to_x[], and also 
    returns the alignment score as the return value.
*/

LPOScore_T align_lpo_po (LPOSequence_T *lposeq_x,
			 LPOSequence_T *lposeq_y,
			 ResidueScoreMatrix_T *m,
			 LPOLetterRef_T **x_to_y,
			 LPOLetterRef_T **y_to_x,
			 LPOScore_T (*scoring_function)
			 (int, int, LPOLetter_T *, LPOLetter_T *, ResidueScoreMatrix_T *),
			 int use_global_alignment)
{
  LPOLetter_T *seq_x = lposeq_x->letter;
  LPOLetter_T *seq_y = lposeq_y->letter;
  
  int len_x, len_y;
  int n_edges_x, n_edges_y;
  int *node_type_x, *node_type_y;
  int *refs_from_right_x, *refs_from_right_y;
  int max_rows_alloced_x, max_rows_alloced_y, n_score_rows_alloced = 0;
  
  int i, j, xcount, ycount, prev_gap;
  int best_x = -1, best_y = -1;
  LPOScore_T min_score = -999999, best_score = -999999;
  int possible_end_square;
  LPOLetterLink_T **x_left = NULL, **y_left = NULL, *xl, *yl;
  DPMove_T **move = NULL, *my_move;
  
  DPScore_T *curr_score = NULL, *prev_score = NULL, *init_col_score = NULL, *my_score;
  DPScore_T **score_rows = NULL;

  int max_gap_length;
  LPOScore_T *gap_penalty_x, *gap_penalty_y;
  int *next_gap_array, *next_perp_gap_array;
  
  LPOScore_T try_score, insert_x_score, insert_y_score, match_score;
  int insert_x_x, insert_x_gap;
  int insert_y_y, insert_y_gap;
  int match_x, match_y;
  
  get_lpo_stats (lposeq_x, &len_x, &n_edges_x, &node_type_x, &refs_from_right_x, &max_rows_alloced_x, &x_left);
  get_lpo_stats (lposeq_y, &len_y, &n_edges_y, &node_type_y, &refs_from_right_y, &max_rows_alloced_y, &y_left);

  /*
    fprintf (stdout, "sequence x:  %ld nodes, %ld edges, %ld rows at most --> %ld mem\n", len_x, n_edges_x, max_rows_alloced_x, max_rows_alloced_x * len_y);
    fprintf (stdout, "sequence y:  %ld nodes, %ld edges, %ld rows at most --> %ld mem\n", len_y, n_edges_y, max_rows_alloced_y, max_rows_alloced_y * len_x);
  */
  
  /* INITIALIZE GAP PENALTIES: */
  max_gap_length = m->max_gap_length;
  gap_penalty_x = m->gap_penalty_x;
  gap_penalty_y = m->gap_penalty_y;
  CALLOC (next_gap_array, max_gap_length + 2, int);
  CALLOC (next_perp_gap_array, max_gap_length + 2, int);

  for (i=0; i<max_gap_length+1; i++) {
    /* GAP LENGTH EXTENSION RULE: */
    /* 0->1, 1->2, 2->3, ..., M-1->M; but M->M. */
    next_gap_array[i] = (i<max_gap_length) ? i+1 : i;
    /* PERPENDICULAR GAP (i.e. X FOR A GROWING Y-GAP) IS KEPT AT 0 IF DOUBLE-GAP-SCORING (old scoring) IS USED. */
    next_perp_gap_array[i] = (DOUBLE_GAP_SCORING ? 0 : next_gap_array[i]);
  }
  
  /* GAP LENGTH = M+1 IS USED FOR INITIAL STATE. */
  /* THIS MUST BE TREATED DIFFERENTLY FOR GLOBAL v. LOCAL ALIGNMENT: */
  if (0 == use_global_alignment) {   /* FREE EXTENSION OF INITIAL GAP (FOR LOCAL ALIGNMENT) */
    gap_penalty_x[max_gap_length+1] = gap_penalty_y[max_gap_length+1] = 0;
    next_gap_array[max_gap_length+1] = next_perp_gap_array[max_gap_length+1] = max_gap_length+1;
  }
  else {   /* TREAT INITIAL GAP LIKE ANY OTHER (FOR GLOBAL ALIGNMENT) */
    gap_penalty_x[max_gap_length+1] = gap_penalty_x[0];
    gap_penalty_y[max_gap_length+1] = gap_penalty_y[0];
    next_gap_array[max_gap_length+1] = next_gap_array[0];
    next_perp_gap_array[max_gap_length+1] = next_perp_gap_array[0];
  }
  
  
  /* ALLOCATE MEMORY FOR 'MOVE' AND 'SCORE' MATRICES: */
  
  CALLOC (move, len_y, DPMove_T *);
  for (i=0; i<len_y; i++) {
    CALLOC (move[i], len_x, DPMove_T);
  }

  CALLOC (init_col_score, len_y+1, DPScore_T);
  init_col_score = &(init_col_score[1]);
  
  CALLOC (score_rows, len_y+1, DPScore_T *);
  score_rows = &(score_rows[1]);
  CALLOC (score_rows[-1], len_x+1, DPScore_T);
  score_rows[-1] = &(score_rows[-1][1]);
  curr_score = score_rows[-1];


  /* FILL INITIAL ROW (-1). */
  /* GAP LENGTH = M+1 IS USED FOR INITIAL STATE. */
  
  curr_score[-1].score = 0;
  curr_score[-1].gap_x = curr_score[-1].gap_y = max_gap_length+1;
  
  for (i=0; i<len_x; i++) {
    curr_score[i].score = min_score;
    for (xcount = 1, xl = x_left[i]; xl != NULL; xcount++, xl = xl->more) {
      prev_gap = curr_score[xl->ipos].gap_x;
      try_score = curr_score[xl->ipos].score + xl->score - gap_penalty_x[prev_gap];
      if (try_score > curr_score[i].score) {
	curr_score[i].score = try_score;
	curr_score[i].gap_x = next_gap_array[prev_gap];
	curr_score[i].gap_y = next_perp_gap_array[prev_gap];
      }
    }
  }
  
  /* FILL INITIAL COLUMN (-1). */
  
  init_col_score[-1] = curr_score[-1];
  for (i=0; i<len_y; i++) {
    init_col_score[i].score = min_score;
    for (ycount = 1, yl = y_left[i]; yl != NULL; ycount++, yl = yl->more) {
      prev_gap = init_col_score[yl->ipos].gap_y;
      try_score = init_col_score[yl->ipos].score + yl->score - gap_penalty_y[prev_gap];
      if (try_score > init_col_score[i].score) {
	init_col_score[i].score = try_score;
	init_col_score[i].gap_x = next_perp_gap_array[prev_gap];
	init_col_score[i].gap_y = next_gap_array[prev_gap];
      }
    }
  }

  
  /** MAIN DYNAMIC PROGRAMMING LOOP **/


  /* OUTER LOOP (i-th position in LPO y): */
  for (i=0; i<len_y; i++) {
    
    /* ALLOCATE MEMORY FOR 'SCORE' ROW i: */
    CALLOC (score_rows[i], len_x+1, DPScore_T);
    score_rows[i] = &(score_rows[i][1]);
    n_score_rows_alloced++;
        
    curr_score = score_rows[i];
    curr_score[-1] = init_col_score[i];
          
    /* INNER LOOP (j-th position in LPO x): */
    for (j=0; j<len_x; j++) {

      match_score = (use_global_alignment) ? min_score : 0;
      match_x = match_y = 0;
      
      insert_x_score = insert_y_score = min_score;
      insert_x_x = insert_y_y = 0;
      insert_x_gap = insert_y_gap = 0;
      
      /* THIS SQUARE CAN END THE ALIGNMENT IF WE'RE USING LOCAL ALIGNMENT, */
      /* OR IF BOTH THE X- AND Y-NODES CONTAIN THE END OF A SEQUENCE. */
      possible_end_square = ((0 == use_global_alignment) || ((node_type_x[j] & LPO_FINAL_NODE) && (node_type_y[i] & LPO_FINAL_NODE)));
      
      /* LOOP OVER y-predecessors: */
      for (ycount = 1, yl = y_left[i]; yl != NULL; ycount++, yl = yl->more) {
	
	prev_score = score_rows[yl->ipos];
	
	/* IMPROVE Y-INSERTION?: trace back to (i'=yl->ipos, j) */
	prev_gap = prev_score[j].gap_y;
	try_score = prev_score[j].score + yl->score - gap_penalty_y[prev_gap];
	if (try_score > insert_y_score) {
	  insert_y_score = try_score;
	  insert_y_y = ycount;
	  insert_y_gap = prev_gap;
	}
	
	/* LOOP OVER x-predecessors (INSIDE y-predecessor LOOP): */
	for (xcount = 1, xl = x_left[j]; xl != NULL; xcount++, xl = xl->more) {
	  
	  /* IMPROVE XY-MATCH?: trace back to (i'=yl->ipos, j'=xl->ipos) */
	  try_score = prev_score[xl->ipos].score + xl->score + yl->score;
	  if (try_score > match_score) {
	    match_score = try_score;
	    match_x = xcount;
	    match_y = ycount;
	  }
	}
      }
      
      /* LOOP OVER x-predecessors (OUTSIDE y-predecessor LOOP): */
      for (xcount = 1, xl = x_left[j]; xl != NULL; xcount++, xl = xl->more) {

	/* IMPROVE X-INSERTION?: trace back to (i, j'=xl->ipos) */
	prev_gap = curr_score[xl->ipos].gap_x;
	try_score = curr_score[xl->ipos].score + xl->score - gap_penalty_x[prev_gap];
	if (try_score > insert_x_score) {
	  insert_x_score = try_score;
	  insert_x_x = xcount;
	  insert_x_gap = prev_gap;
	}
      }
      
      /* USE CUSTOM OR DEFAULT SCORING FUNCTION: */
      if (scoring_function != NULL) {
	match_score += scoring_function (j, i, seq_x, seq_y, m);
      }
      else {
	match_score += m->score[seq_x[i].letter][seq_y[j].letter];
      }
      
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

      /* RECORD BEST ALIGNMENT END FOR TRACEBACK: */
      if (possible_end_square && my_score->score >= best_score) {
	/* BREAK TIES BY CHOOSING MINIMUM (x,y): */
	if (my_score->score > best_score || (j == best_x && i < best_y) || j < best_x) {
	  best_score = my_score->score;
	  best_x = j;
	  best_y = i;
	}
      }
    }

    /* UPDATE # OF REFS TO 'SCORE' ROWS; FREE MEMORY WHEN POSSIBLE: */
    for (yl = y_left[i]; yl != NULL; yl = yl->more) if ((j = yl->ipos) >= 0) {
      if ((--refs_from_right_y[j]) == 0) {
	score_rows[j] = &(score_rows[j][-1]);
	FREE (score_rows[j]);
	n_score_rows_alloced--;
      }
    }
    if (refs_from_right_y[i] == 0) {
      score_rows[i] = &(score_rows[i][-1]);
      FREE (score_rows[i]);
      n_score_rows_alloced--;
    }
  }
  
  IF_GUARD(best_x>=len_x || best_y>=len_y,1.1,(ERRTXT,"Bounds exceeded!\nbest_x,best_y:%d,%d\tlen:%d,%d\n",best_x,best_y,len_x,len_y),CRASH);
  
  /**/
    fprintf (stderr, "aligned (%d nodes, %ld edges) to (%d nodes, %ld edges): ", len_x, n_edges_x, len_y, n_edges_y);
    fprintf (stderr, "best %s score = %d @ (%d %d)\n", (use_global_alignment ? "global" : "local"), best_score, best_x, best_y);
    /**/
    
  /* DYNAMIC PROGRAMING MATRIX COMPLETE, NOW TRACE BACK FROM best_x, best_y */
  trace_back_lpo_alignment (len_x, len_y, move, x_left, y_left,
			    best_x, best_y,
			    x_to_y, y_to_x);


  /* CLEAN UP AND RETURN: */
  
  FREE (node_type_x);
  FREE (node_type_y);
  
  FREE (refs_from_right_x);
  FREE (refs_from_right_y);

  FREE (next_gap_array);
  FREE (next_perp_gap_array);
  
  score_rows[-1] = &(score_rows[-1][-1]);
  FREE (score_rows[-1]);
  score_rows = &(score_rows[-1]);
  FREE (score_rows);
  
  init_col_score = &(init_col_score[-1]);
  FREE (init_col_score);
    
  for (i=0; i<len_x; i++) {
    if (x_left[i] != &seq_x[i].left) {
      FREE (x_left[i]);
    }
  }
  FREE (x_left);
  
  for (i=0; i<len_y; i++) {
    if (y_left[i] != &seq_y[i].left) {
      FREE (y_left[i]);
    }
  }
  FREE (y_left);
  
  for (i=0; i<len_y; i++) {
    FREE (move[i]);
  }
  FREE (move);
  
  return best_score;
}

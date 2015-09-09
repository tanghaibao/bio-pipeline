
#include "default.h"
#include "seq_util.h"



/** randomizes seq[] by shuffling, and places the result in randseq[];
 if randseq[] and seq[] are distinct, seq[] is left unchanged */
void shuffle_seq(int len,
		char seq[],
		char randseq[])
{
  int i,j;
  char c;

  for (i=0;i<len;i++)  {
    j=rand()%len; /* CHOOSE RANDOM POSITION j */
    c=seq[i]; /* SWAP POSITIONS i AND j */
    randseq[i]=seq[j];
    randseq[j]=c;
  }

  return;
}






/**@memo TRANSLATE FROM ASCII LETTERS TO COMPACTED NUMBERICAL INDEX: 
    index_symbols(seq[i].length,seq[i].sequence,seq[i].sequence,
		  m->nsymbol,m->symbol);
*/
/** converts characters in seq[] to the INDEX of the matching character in
 symbols[], and returns the result in out[] */
void index_symbols(int nseq,char seq[],char out[],
		   int nsymbs,char symbols[])
{
  int i,j,k;
  LOOP (i,nseq) {
    k=nsymbs-1; /* DEFAULT: UNMATCHABLE SYMBOL */
    LOOP (j,nsymbs) {  /* FIND MATCHING SYMBOL */
      if (symbols[j]==seq[i]) {  /* FOUND IT! */
	k=j;
	break;
      }
    }
    out[i]=k; /* SAVE THE TRANSLATED CODE */
  }
  return;
}






int *Score_matrix_row=NULL;

int best_match_qsort_cmp(const void *void_a,const void *void_b)
{
  int *a=(int *)void_a,*b=(int *)void_b;

  if (Score_matrix_row[*a]>Score_matrix_row[*b])
    return -1;
  else if (Score_matrix_row[*a]<Score_matrix_row[*b])
    return 1;
  else /* EQUAL */
    return 0;
}




#ifdef SOURCE_EXCLUDED
char DNA_symbols[1024];
float DNA_rescale_score;
#endif

/** reads an alignment scoring matrix in the pam format */
int read_score_matrix(char filename[],ResidueScoreMatrix_T *m)
{
  int i,j,k,nsymb=0,found_symbol_line=0,isymb;
  char line[1024],dna_codes[256];
  FILE *ifile;

   /* GAP PENALTY DEFAULTS */
  m->gap_penalty_set[0][0]=m->gap_penalty_set[1][0]=12; /*SAVE PENALTIES*/
  m->gap_penalty_set[0][1]=m->gap_penalty_set[1][1]=2;
  m->gap_penalty_set[0][2]=m->gap_penalty_set[1][2]=0;
  m->trunc_gap_length = TRUNCATE_GAP_LENGTH;
  m->decay_gap_length = DECAY_GAP_LENGTH;
  
  ifile=fopen(filename,"r");
  if (!ifile) {
    WARN_MSG(USERR,(ERRTXT,"Can't open alignment matrix from %s\n",filename),"$Revision: 1.2.2.2 $");
    return -2; /* FAILED TO FIND FILE TO READ */
  }

  while (fgets(line,1023,ifile)) {
    if ('#'==line[0] || '\n'==line[0]) /* SKIP COMMENT OR BLANK LINES */
      continue;
    
    else if (1==sscanf(line,"GAP-TRUNCATION-LENGTH=%d",&i)) {
      m->trunc_gap_length = i;
    }
    
    else if (1==sscanf(line,"GAP-DECAY-LENGTH=%d",&i)) {
      m->decay_gap_length = i;
    }
    
    else if (3==sscanf(line,"GAP-PENALTIES=%d %d %d",&i,&j,&k)) {
      m->gap_penalty_set[0][0]=m->gap_penalty_set[1][0]=i; /*SAVE PENALTIES*/
      m->gap_penalty_set[0][1]=m->gap_penalty_set[1][1]=j;
      m->gap_penalty_set[0][2]=m->gap_penalty_set[1][2]=k;
    }

    else if (3==sscanf(line,"GAP-PENALTIES-X=%d %d %d",&i,&j,&k)) {
      m->gap_penalty_set[1][0]=i; /*SAVE PENALTIES ONLY FOR X DIRECTION*/
      m->gap_penalty_set[1][1]=j;
      m->gap_penalty_set[1][2]=k;
    }

#ifdef SOURCE_EXCLUDED
    else if (1==sscanf(line,"DNACODES=%99s",dna_codes)) { /* READ DNACODES*/
      strcpy(DNA_symbols,dna_codes);/*SYMBOLS COUNTED AS DNA FOR AUTORECOG*/
    }

    else if (1==sscanf(line,"DNASCALE=%f",&DNA_rescale_score))
      continue;
#endif
    
    else if (!found_symbol_line) { /* READ THIS LINE AS LIST OF SEQ SYMBOLS*/
      for (i=0;'\0'!=line[i];i++)
	if (!isspace(line[i])) /* IGNORE WHITESPACE */
	  m->symbol[nsymb++]=line[i]; /* SAVE TO LIST OF SYMBOLS */
      found_symbol_line=1; /* SET FLAG SO WE NOW READ MATRIX SCORE VALUES */
    }
    
    else { /* READ SCORING MATRIX LINES */
      found_symbol_line=0; /* DEFAULT: FAILED TO FIND MATCHING SYMBOL IN LIST*/
      LOOP (isymb,nsymb) /* FIND MATCH TO THIS SYMBOL */
	if (m->symbol[isymb]==line[0]) {
	  found_symbol_line=1; /* SIGNAL THAT WE SUCCESFULLY FOUND MATCH */
	  j=1; /* SKIP FIRST CHARACTER: OUR SEQUENCE SYMBOL */
	  LOOPF (i,nsymb) { /* READ ALL THE SCORE VALUES ON THIS LINE */
	    if (1==sscanf(line+j,"%d%n",&(m->score[isymb][i]),&k))
	      j+=k; /* ADVANCE THE READING POSITION */
	    else { /* MISSING SCORE DATA: ERROR! */
	      IF_GUARD(1,5.23,(ERRTXT,"Missing score value for pair %c:%c",
			      m->symbol[isymb],m->symbol[i]),TRAP)
		;
	      fclose(ifile); /* CLOSE OUR STREAM */
	      return -1;
	    }
	  }
	  break;
	}
      IF_GUARD(!found_symbol_line,1.5,(ERRTXT,"Missing or unknown sequence symbol: %c",line[0]),TRAP) { /* ERROR: AN INVALID SYMBOL, NOT IN LIST */
	fclose(ifile); /* CLOSE OUR STREAM */
	return -1;
      }
    }
  }
  fclose(ifile);

  /* CONSTRUCT GAP PENALTY ARRAYS FROM GAP PARAMETERS: */
  m->max_gap_length = m->trunc_gap_length + m->decay_gap_length;
  CALLOC (m->gap_penalty_x, m->max_gap_length+2, LPOScore_T);
  CALLOC (m->gap_penalty_y, m->max_gap_length+2, LPOScore_T);

  /*** GAP OPENING PENALTY @ L=0->1 */
  m->gap_penalty_x[0] = m->gap_penalty_set[0][0];
  m->gap_penalty_y[0] = m->gap_penalty_set[1][0];
  
  /*** 1st AFFINE EXTENSION PENALTY (A1) @ L=1->2,2->3,...T-1->T */
  for (i=1;i<m->trunc_gap_length;i++) {
    m->gap_penalty_x[i] = m->gap_penalty_set[0][1];
    m->gap_penalty_y[i] = m->gap_penalty_set[1][1];
  }
  
  /*** DECAYING EXTENSION PENALTY (A1-->A2; skipped if D=0) @ L=T->T+1,...T+D-1->T+D */
  for (i=0;i<m->decay_gap_length;i++) {
    double dec_x = (m->gap_penalty_set[0][1] - m->gap_penalty_set[0][2]) / ((double)(m->decay_gap_length + 1));
    double dec_y = (m->gap_penalty_set[1][1] - m->gap_penalty_set[1][2]) / ((double)(m->decay_gap_length + 1));
    m->gap_penalty_x[i+m->trunc_gap_length] = m->gap_penalty_set[0][1] - (i+1) * dec_x;
    m->gap_penalty_y[i+m->trunc_gap_length] = m->gap_penalty_set[1][1] - (i+1) * dec_y;
  }
  
  /*** 2nd AFFINE EXTENSION PENALTY (A2) @ L>=T+D */
  m->gap_penalty_x[m->max_gap_length] = m->gap_penalty_set[0][2];
  m->gap_penalty_y[m->max_gap_length] = m->gap_penalty_set[1][2];
  
  m->gap_penalty_x[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */
  m->gap_penalty_y[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */
  
  
  LOOPF (i,nsymb) {
    Score_matrix_row= m->score[i]; /* ROW TO USE FOR SORTING best_match */
    LOOP (j,nsymb)
      m->best_match[i][j] = j;
    qsort(m->best_match[i],nsymb,sizeof(m->best_match[0][0]),
	  best_match_qsort_cmp);
#ifdef SOURCE_EXCLUDED
    printf("%c SORT",m->symbol[i]); /* TEST: PRINT OUT SORTED TABLE */
    LOOPF (j,nsymb)
      printf("\t%c:%d",m->symbol[m->best_match[i][j]],
	     m->score[i][m->best_match[i][j]]);
    printf("\n");
#endif
  }

  m->symbol[nsymb]='\0'; /* TERMINATE THE SYMBOL STRING */
  m->nsymbol=nsymb;
  return nsymb;
}



/** prints a scoring matrix, only including those symbols in subset[] */
void print_score_matrix(FILE *ifile,ResidueScoreMatrix_T *m,char subset[])
{
  int i,i_m,j,j_m,nsubset;

  nsubset=strlen(subset);

  printf(" ");
  LOOPF (i,nsubset)
    printf("  %c",subset[i]);
  printf("\n");

  LOOPF (i,nsubset) {
    LOOP (i_m,m->nsymbol)
      if (m->symbol[i_m]==subset[i])
	break;
    printf("%c",subset[i]);
    LOOPF (j,nsubset) {
      LOOP (j_m,m->nsymbol)
	if (m->symbol[j_m]==subset[j])
	  break;
      printf("%3d",m->score[i_m][j_m]);
    }
    printf("\n");
  }
  return;
}



/** restricts seq[] to the set of allowed characters given by symbol[];
 other characters will be replaced by the default symbol[0] */
int limit_residues(char seq[],char symbol[])
{
  int i,len,nreplace=0;

  len=strlen(seq);
  for (i=strspn(seq,symbol);i<len;i=i+1+strspn(seq+i+1,symbol)) {
    seq[i]=symbol[0]; /* FORCE IT TO BE A VALID SYMBOL */
    nreplace++; /* COUNT TOTAL REPLACEMENTS */
  }
  return nreplace;
}





/** RETURNS THE COMPLEMENTARY BASE *IN LOWER CASE* */
char complementary_base(char base) 
{
  switch (base) {
  case 'A': case 'a':                     return 't';
  case 'T': case 't': case 'U': case 'u': return 'a';
  case 'G': case 'g':                     return 'c';
  case 'C': case 'c':                     return 'g';
  default: return base;
  }
}

/** REVERSE COMPLEMENTS seq[] IN PLACE, AND RETURNS POINTER TO seq  
---------------------------------------------------------------
-----------------------------------------------------------
*/
char *reverse_complement(char seq[])
{
  int i,j;
  char c;
  for ((i=0),(j=strlen(seq)-1);i<=j;i++,j--) {/* SWAP FROM ENDS TO CENTER*/
    c=complementary_base(seq[i]); /* SWAP ENDS AND COMPLEMENT */
    seq[i]=complementary_base(seq[j]);
    seq[j]=c;
  }
  return seq;
}



#ifndef SEQ_UTIL_HEADER_INCLUDED
#define SEQ_UTIL_HEADER_INCLUDED

/* NB: YOU *MUST* TYPEDEF Sequence_T ELSEWHERE, IN ORDER FOR THE PROTOTYPES
   BELOW TO WORK... */

/* USE_PROJECT_HEADER IS A GENERAL WAY TO SNEAK IN CUSTOMIZATION
   BEFORE PROCESSING DEFINITIONS IN THIS FILE */
#ifdef USE_PROJECT_HEADER
#include "project.h"
#endif

#ifndef RESIDUE_SCORE_DEFINED
typedef int ResidueScore_T; /* DEFAULT: USE int FOR SCORING */
#endif




#define FASTA_NAME_MAX   4096  /*define the max length of the name */
#define SEQ_LENGTH_MAX  32768 /* define the maximum length of each seq */

enum {
  dont_switch_case,
  switch_case_to_lower,
  switch_case_to_upper
};



#define MATRIX_SYMBOL_MAX 128
typedef struct {
  int nsymbol;
  char symbol[MATRIX_SYMBOL_MAX];
  ResidueScore_T score[MATRIX_SYMBOL_MAX][MATRIX_SYMBOL_MAX];
  int best_match[MATRIX_SYMBOL_MAX][MATRIX_SYMBOL_MAX];
  
  ResidueScore_T gap_penalty_set[2][3]; 
  int trunc_gap_length;
  int decay_gap_length;
  
  ResidueScore_T *gap_penalty_x, *gap_penalty_y;
  int max_gap_length;
  
  int nfreq; /* STORE FREQUENCIES OF AMINO ACIDS FOR BALANCING MATRIX...*/
  char freq_symbol[MATRIX_SYMBOL_MAX];
  float freq[MATRIX_SYMBOL_MAX];
} ResidueScoreMatrix_T;

/****************************************************** seq_util.c */
void shuffle_seq(int len,
		char seq[],
		char randseq[]);

void index_symbols(int nseq,char seq[],char out[],
		   int nsymbs,char symbols[]);

int read_score_matrix(char filename[],ResidueScoreMatrix_T *m);

void print_score_matrix(FILE *ifile,ResidueScoreMatrix_T *m,char subset[]);

int limit_residues(char seq[],char symbol[]);

int create_seq(int nseq,Sequence_T **seq,char seq_name[],char seq_title[],char tmp_seq[],int do_switch_case);

char *reverse_complement(char seq[]);



/******************************************************* fasta_format.c */
int read_fasta(FILE *seq_file,Sequence_T **seq,
	       int do_switch_case,char **comment);

void write_fasta(FILE *ifile,char name[],char title[],char seq[]);

#endif

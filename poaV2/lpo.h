

#ifndef LPO_HEADER_INCLUDED
#define LPO_HEADER_INCLUDED

#include <default.h>
#include <poa.h>
#include <seq_util.h>

/*********************************************************** lpo.c */
void lpo_init(LPOSequence_T *seq);

void initialize_seqs_as_lpo(int nseq, Sequence_T seq[],ResidueScoreMatrix_T *m);
void lpo_index_symbols(LPOSequence_T *lpo,ResidueScoreMatrix_T *m);

/** reindex source sequences in `seq' so that the sequence in position i
    ends up in position perm[i].  `perm' must be a permutation of the
    integers [0,seq->nsource_seq-1]. */
void reindex_lpo_source_seqs (LPOSequence_T *seq, int *perm);

int save_lpo_source(LPOSequence_T *seq,
		     char name[],
		     char title[],
		     int length,
		     int weight,
		     int bundle_id,
		    int ndata,
		    LPONumericData_T data[]);

int *save_lpo_source_list(LPOSequence_T *seq,
			  int nsource_seq,
			  LPOSourceInfo_T source_seq[]);

LPOLetterLink_T *add_lpo_link(LPOLetterLink_T *list,LPOLetterRef_T ipos);

void add_lpo_sources(LPOLetterSource_T *new_s,LPOLetterSource_T *old_s,
		     int iseq_new[]);

void crosslink_rings(LPOLetterRef_T a,LPOLetterRef_T b,LPOLetter_T seq[]);

LPOSequence_T *copy_fuse_lpo(LPOSequence_T *holder_x,
			     LPOSequence_T *holder_y,
			     LPOLetterRef_T x_to_y[],
			     LPOLetterRef_T y_to_x[]);

LPOSequence_T *copy_lpo(LPOSequence_T *holder_x);

LPOSequence_T *fuse_lpo(LPOSequence_T *holder_x,
			LPOSequence_T *holder_y,
			LPOLetterRef_T x_to_y[],
			LPOLetterRef_T y_to_x[]);

void free_lpo_letters(int nletter,LPOLetter_T *letter,int please_free_block);

void free_lpo_sequence(LPOSequence_T *seq,int please_free_holder);

int add_path_sequence(int path_length,
		      LPOLetterRef_T path[],
		      LPOSequence_T *seq,
		      char name[],
		      char title[]);


/************************************************** FROM align_lpo.c */
int align_lpo(LPOSequence_T *lposeq_x,
	      LPOSequence_T *lposeq_y,
	      ResidueScoreMatrix_T *m,
	      LPOLetterRef_T **x_to_y,
	      LPOLetterRef_T **y_to_x,
	      int use_global_alignment);



/************************************************** FROM align_lpo_po.c */
LPOScore_T align_lpo_po(LPOSequence_T *lposeq_x,
			LPOSequence_T *lposeq_y,
			ResidueScoreMatrix_T *m,
			LPOLetterRef_T **x_to_y,
			LPOLetterRef_T **y_to_x,
			LPOScore_T (*scoring_function)
			 (int,int,LPOLetter_T [],LPOLetter_T [],
			  ResidueScoreMatrix_T *),
			int use_global_alignment);


/************************************************** FROM buildup_lpo.c */
LPOSequence_T *buildup_lpo(LPOSequence_T *new_seq,
			   int nseq,LPOSequence_T seq[],
			   ResidueScoreMatrix_T *score_matrix,
			   int use_aggressive_fusion,
			   int use_global_alignment);

LPOSequence_T *buildup_clipped_lpo(LPOSequence_T *new_seq,
				  int nseq,LPOSequence_T seq[],
				  ResidueScoreMatrix_T *score_matrix,
				   int use_global_alignment);

LPOSequence_T *buildup_progressive_lpo(int nseq, LPOSequence_T **seqs,
				       ResidueScoreMatrix_T *score_matrix,
				       int use_aggressive_fusion,
                                       int do_progressive,
				       char score_file[], 
				       LPOScore_T (*scoring_function)
				       (int,int,LPOLetter_T [],LPOLetter_T [],
					ResidueScoreMatrix_T *),
                                       int use_global_alignment,
				       int preserve_sequence_order);
				       
LPOSequence_T *buildup_pairwise_lpo(LPOSequence_T seq1[],LPOSequence_T seq2[],
				    ResidueScoreMatrix_T *score_matrix,
				    int use_aggressive_fusion,
				    LPOScore_T (*scoring_function)
				    (int,int,LPOLetter_T [],LPOLetter_T [],
				     ResidueScoreMatrix_T *),
                                    int use_global_alignment);
				    
/**************************************************** lpo_format.c */
void write_lpo(FILE *ifile,LPOSequence_T *seq,
	       ResidueScoreMatrix_T *score_matrix);

LPOSequence_T *read_lpo(FILE *ifile);
LPOSequence_T *read_lpo_select(FILE *ifile,FILE *select_ifile,
			       int keep_all_links,int remove_listed_sequences);

void write_lpo_as_fasta(FILE *ifile,LPOSequence_T *seq,
			int nsymbol,char symbol[]);

void write_lpo_bundle_as_fasta(FILE *ifile,LPOSequence_T *seq,
			       int nsymbol,char symbol[],int ibundle);
void export_clustal_seqal(FILE *ifile,
			  LPOSequence_T *seq,
			  int nsymbol,char symbol[]);


/****************************************************** heaviest_bundle.c */
void generate_lpo_bundles(LPOSequence_T *seq,float minimum_fraction);


/****************************************************** make_frame.c */
LPOSequence_T *build_3_frames(char dna_seq[],char name[],char title[],
			      LPOScore_T frameshift_score,
			      ResidueScoreMatrix_T *m);
LPOSequence_T *map_protein_to_dna(char dna_name[],
				  LPOSequence_T *lpo_dna,
				  int nseq,
				  Sequence_T seq[],
				  ResidueScoreMatrix_T *score_matrix,
				  int use_aggressive_fusion);


/****************************************************** remove_bundle.c */
int remove_bundle(LPOSequence_T *seq,int ibundle,int delete_all_others);



/******************************************************* numeric_data.c */
LPONumericData_T *new_numeric_data(LPOSourceInfo_T *source_seq,
				   char name[],
				   char title[],
				   double initial_value);

LPONumericData_T *find_numeric_data(LPOSourceInfo_T *source_seq,
				   char name[]);

void free_lpo_numeric_data(int ndata,LPONumericData_T *data,
			   int please_free_block);

void new_numeric_data_sets(LPOSourceInfo_T *source_seq,
			   int nset,char *set_names[],
			   char source_name_fmt[],
			   char target_name_fmt[],
			   char title_fmt[]);
void read_numeric_data(int nsource_seq,
		       LPOSourceInfo_T source_seq[],
		       FILE *ifile);
LPONumericData_T *cp_numeric_data(LPOSourceInfo_T *source_seq,
				  LPONumericData_T *data);

/******************************************************* balance_matrix.c */
int read_aa_frequencies(char filename[],ResidueScoreMatrix_T *score_matrix)
     ;
void balance_matrix_score(int nletter,LPOLetter_T letter[],
			  ResidueScoreMatrix_T *score_matrix);


#endif

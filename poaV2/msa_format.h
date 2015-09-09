

/****************/
/* msa_format.h */
/****************/

/* ---
   Functions for reading CLUSTAL- and FASTA-PIR-formatted files
   into the LPOSequence_T data structure, and for determining file
   type from the initial line(s) of a file.
   --- */


#ifndef MSA_FORMAT_HEADER_INCLUDED
#define MSA_FORMAT_HEADER_INCLUDED

#include "default.h"
#include "seq_util.h"
#include "lpo.h"

/** types of MSA files supported */
typedef enum {
  UNKNOWN_MSA, CLUSTAL_MSA, PIR_MSA, PO_MSA
}
msa_file_format;


/** Reads an MSA from a file (cf. read_msa_select).
 */
LPOSequence_T *read_msa (FILE *ifile, msa_file_format format,
			 int do_switch_case, ResidueScoreMatrix_T *score_matrix);


/** Reads an MSA from a file.  If `format' is UNKNOWN_MSA, the file
    format is determined from the first line(s) of the file.
    Uses `select_ifile', if non-NULL, to filter the sequence set.
    Uppercases or lowercases sequence characters according to `do_switch_case'.
    Indexes LPO symbols using `score_matrix'.
*/
LPOSequence_T *read_msa_select (FILE *ifile, msa_file_format format,
				FILE *select_ifile, int keep_all_links, int remove_listed_sequences,
				int do_switch_case, ResidueScoreMatrix_T *score_matrix);


/** Creates an LPO from an RC-MSA alignment matrix.
 */
LPOSequence_T *lpo_from_aln_mat (int n_seqs, char **seq_names, char **seq_titles, char **aln_mat, int *aln_lengths,
				 int do_switch_case, ResidueScoreMatrix_T *score_matrix);



#endif  /* MSA_FORMAT_HEADER_INCLUDED */

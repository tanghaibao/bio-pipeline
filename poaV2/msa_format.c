

/****************/
/* msa_format.c */
/****************/

/* ---
   Implementation of "msa_format.h":
   
   Functions for reading CLUSTAL- and FASTA-PIR-formatted files
   into the LPOSequence_T data structure, and for determining file
   type from the initial line(s) of a file.
   --- */



#include "msa_format.h"


/** is `ch' an allowed residue? (a-z OR A-Z OR ? OR [ OR ]) */
static int is_residue_char (char ch);

/** is `ch' an allowed gap? (. OR -) */
static int is_gap_char (char ch);

/** could `ch' be the first character of a sequence name? (NOT # AND NOT * AND NOT whitespace) */
static int is_name_first_char (char ch);


/** decide which sequences in an RC-MSA alignment to keep based on a filter file `select_ifile'.
    remove discarded sequences and reindex the remaining ones.
    returns the NEW number of sequences. */
static int filter_sequence_set (int n_seqs, FILE *select_ifile, int remove_listed_sequences,
				char **seq_names, char **seq_titles, char **aln_mat, int *aln_lengths);


/** remove sequence bad_seq_id from each letter containing it.
    ASSUMES THAT NO LETTER CONTAINS _ONLY_ bad_seq_id... e.g. bad_seq_id IS A CONSENSUS.
    either keep or don't keep links present only in the removed sequence. */ 
static void strip_seq_from_lpo (LPOSequence_T *lposeq, int bad_seq_id, int keep_all_links);


static LPOSequence_T *read_clustal (FILE *ifile, const char *first_line,
				    FILE *select_ifile, int remove_listed_sequences,
				    int do_switch_case, ResidueScoreMatrix_T *score_matrix);

static LPOSequence_T *read_pir (FILE *ifile, const char *first_line,
				FILE *select_ifile, int remove_listed_sequences,
				int do_switch_case, ResidueScoreMatrix_T *score_matrix);



/** Reads an MSA from a file.  If `format' is UNKNOWN_MSA, the file
    format is determined from the first line(s) of the file.
    Uses `select_ifile', if non-NULL, to filter the sequence set.
    Uppercases or lowercases sequence characters according to `do_switch_case'.
    Indexes LPO symbols using `score_matrix'.
*/
LPOSequence_T *read_msa_select (FILE *ifile, msa_file_format format,
				FILE *select_ifile, int keep_all_links, int remove_listed_sequences,
				int do_switch_case, ResidueScoreMatrix_T *score_matrix)
{
  static char line[512]="";

  /* USE format TO FIX FILE FORMAT IF POSSIBLE. */
  /* OTHERWISE, PO FILES START WITH 'VERSION=' (and this line is discarded),
     FASTA-PIR FILES START WITH '>', AND CLUSTAL FILES WITH 'CLUSTAL' OR SIMPLY
     WITH THE FIRST ALIGNMENT LINE.
     LINES STARTING WITH whitespace OR '#' OR '*' ARE IGNORED. */
  
  while (format!=UNKNOWN_MSA || fgets (line, sizeof(line)-1, ifile)) {
    
    if (format==PIR_MSA || line[0] == '>') {
      return read_pir (ifile, line,
		       select_ifile, remove_listed_sequences,
		       do_switch_case, score_matrix);
    }
    else if (format==PO_MSA || 0==strncmp(line,"VERSION=",8)) {
      if (select_ifile != NULL) {
	return read_lpo_select (ifile, select_ifile, keep_all_links, remove_listed_sequences);
      }
      else {
	return read_lpo (ifile);
      }
    }
    else if (format==CLUSTAL_MSA || 0==strncmp(line,"CLUSTAL",7)) {
      return read_clustal (ifile, line,
			   select_ifile, remove_listed_sequences,
			   do_switch_case, score_matrix);
    }
    else if (line[0] == '#' || line[0] == '*' || line[0] == ' ' || line[0] == '\t' || line[0] == '\n' || line[0] == '\r') {
      continue;
    }
    else {
      WARN_MSG(USERR,(ERRTXT, "Unable to determine MSA file type... trying CLUSTAL.\n"),"$Revision: 1.1.2.3 $");
      return read_clustal (ifile, line,
			   select_ifile, remove_listed_sequences,
			   do_switch_case, score_matrix);
    }
  }
  
  WARN_MSG(USERR,(ERRTXT, "No data in MSA file.\n"),"$Revision: 1.1.2.3 $");

  return NULL;
}

/** Reads an MSA from a file (cf. read_msa_select).
 */
LPOSequence_T *read_msa (FILE *ifile, msa_file_format format,
			 int do_switch_case, ResidueScoreMatrix_T *score_matrix)
{
  return read_msa_select (ifile, format,
			  NULL, 0, 0,
			  do_switch_case, score_matrix);
}




/** Reads a CLUSTAL-formatted alignment file.
 */
LPOSequence_T *read_clustal (FILE *fp, const char *first_line,
			     FILE *select_ifile, int remove_listed_sequences,
			     int do_switch_case, ResidueScoreMatrix_T *score_matrix)
{
  int i, j, n_seqs=0, curr_seq=0, expect_repeats=0, expect_header=1, line_num=0;
  char ch;  
  char **seq_names=NULL, **seq_titles=NULL, **aln_mat=NULL;
  int *aln_lengths=NULL;
  char line[512]="", name[512]="", aln[512]="";
  LPOSequence_T *lposeq = NULL;
  
  while ((first_line!=NULL && line_num==0) || fgets (line, sizeof(line)-1, fp)) {
    if (first_line!=NULL && line_num==0) {
      strcpy (line, first_line);
    }
    line_num++;
    
    if (expect_header) {  /* LOOKING FOR 'CLUSTAL' HEADER LINE */
      if (0 == strncmp(line,"CLUSTAL",7)) {  /* FOUND IT; STOP LOOKING */
	expect_header=0;
	continue;
      }
      else if (0 == is_name_first_char(line[0])) {  /* COMMENT LINE; KEEP LOOKING */
	continue;
      }
      else {
	expect_header=0;  /* FOUND SOMETHING ELSE; NO HEADER, SO CONTINUE */
      }
    }
    
    if (0 == is_name_first_char(line[0])) {  /* BLOCK SEPARATOR? */
      curr_seq = 0;
      if (n_seqs>0) {
	expect_repeats=1;  /* YES; NOW EXPECT SAME SEQUENCES OVER AGAIN IN EACH BLOCK */
      }
      continue;
    }
    
    if (sscanf (line, "%s %[^\n\r]", name, aln) < 2) {
      WARN_MSG(USERR,(ERRTXT, "Error: Trouble reading CLUSTAL-formatted file near line %d: \n>>>\n%s<<<\nBailing out.\n",line_num,line),
	       "$Revision: 1.1.2.3 $");
      goto free_memory_and_exit;
    }
    
    if (0 == expect_repeats) {  /* FIRST BLOCK STILL, SO MAKE ROOM FOR NEW SEQ */
      n_seqs++;
      
      REALLOC (seq_names, n_seqs, char *);
      seq_names[n_seqs-1] = strdup(name);
      REALLOC (seq_titles, n_seqs, char *);
      seq_titles[n_seqs-1] = strdup("");
      
      REALLOC (aln_mat, n_seqs, char *);
      aln_mat[n_seqs-1] = strdup(" ");
      REALLOC (aln_lengths, n_seqs, int);
      aln_lengths[n_seqs-1] = 0;
    }
    else if (curr_seq>=n_seqs || strcmp(name,seq_names[curr_seq])) {  /* NAME SHOULD BE A REPEAT */
      WARN_MSG(USERR,(ERRTXT, "Error: Trouble reading CLUSTAL-formatted file at line %d: \n>>>\n%s<<<\nSequence name (%s) does not match expected sequence name (%s).  Bailing out.\n",line_num,line,name,seq_names[curr_seq]),
	       "$Revision: 1.1.2.3 $");
      goto free_memory_and_exit;
    }
    
    REALLOC (aln_mat[curr_seq], aln_lengths[curr_seq] + strlen(aln), char);
    for (i=0; i<strlen(aln); i++) {
      ch = aln[i];
      if (is_residue_char(ch) || is_gap_char(ch)) {
	aln_mat[curr_seq][aln_lengths[curr_seq]++] = aln[i];
      }
    }
    curr_seq++;
  }

  n_seqs = filter_sequence_set (n_seqs, select_ifile, remove_listed_sequences, seq_names, seq_titles, aln_mat, aln_lengths);
  
  lposeq = lpo_from_aln_mat (n_seqs, seq_names, seq_titles, aln_mat, aln_lengths, do_switch_case, score_matrix);
  
 free_memory_and_exit:
  for (i=0; i<n_seqs; i++) {
    FREE (seq_names[i]);
    FREE (seq_titles[i]);
    FREE (aln_mat[i]);
  }
  FREE (seq_names);
  FREE (seq_titles);
  FREE (aln_mat);
  FREE (aln_lengths);
  
  if (n_seqs>0) {
    strcpy (lposeq->name, lposeq->source_seq[0].name);
    FREE (lposeq->title);
    lposeq->title = strdup (lposeq->source_seq[0].title);
  }
  
  return lposeq;
}

/** Reads a FASTA-PIR-formatted alignment file.
 */
LPOSequence_T *read_pir (FILE *fp, const char *first_line,
			 FILE *select_ifile, int remove_listed_sequences,
			 int do_switch_case, ResidueScoreMatrix_T *score_matrix)
{
  int i, j, n_seqs=0, curr_seq=0, line_num=0;
  char ch;  
  char **seq_names=NULL, **seq_titles=NULL, **aln_mat=NULL;
  int *aln_lengths=NULL, *keep_seq=NULL;
  char line[512]="", name[512]="", title[512]="", aln[512]="";
  LPOSequence_T *lposeq = NULL;
  
  while ((first_line!=NULL && line_num==0) || fgets (line, sizeof(line)-1, fp)) {
    if (first_line!=NULL && line_num==0) {
      strcpy (line, first_line);
    }
    line_num++;
    
    if (line[0] == '>') {  /* HEADER LINE FOR NEW SEQUENCE */
      title[0] = '\0';
      if (sscanf (line, ">%s %[^\n\r]", name, title) < 1) {
	WARN_MSG(USERR,(ERRTXT, "Error: Trouble reading PIR-formatted file near line %d (no sequence name?):\n>>>\n%s<<<\nBailing out.\n",line_num,line),
		 "$Revision: 1.1.2.3 $");
	goto free_memory_and_exit;
      }
      
      n_seqs++;
      curr_seq=n_seqs-1;
      REALLOC (seq_names, n_seqs, char *);
      seq_names[n_seqs-1] = strdup(name);
      REALLOC (seq_titles, n_seqs, char *);
      seq_titles[n_seqs-1] = strdup(title);
      REALLOC (aln_mat, n_seqs, char *);
      aln_mat[n_seqs-1] = strdup(" ");
      REALLOC (aln_lengths, n_seqs, int);
      aln_lengths[n_seqs-1] = 0;
    }
    else if (line[0]=='#' || line[0]=='*') {  /* COMMENT LINE */
      continue;
    }
    else {  /* ALIGNMENT ROW FOR CURRENT SEQUENCE */
      sscanf (line, "%[^\n\r]", aln);
      if (n_seqs==0) {
	WARN_MSG(USERR,(ERRTXT, "Error: Trouble reading PIR-formatted file near line %d (no preceding '>seqname' line?):\n>>>\n%s<<<\nBailing out.\n",line_num,line),
		 "$Revision: 1.1.2.3 $");
	goto free_memory_and_exit;
      }
      
      REALLOC (aln_mat[curr_seq], aln_lengths[curr_seq] + strlen(aln), char);
      for (i=0; i<strlen(aln); i++) {
	ch = aln[i];
	if (is_residue_char(ch) || is_gap_char(ch)) {
	  aln_mat[curr_seq][aln_lengths[curr_seq]++] = aln[i];
	}
      }
    }
  }

  n_seqs = filter_sequence_set (n_seqs, select_ifile, remove_listed_sequences, seq_names, seq_titles, aln_mat, aln_lengths);  

  lposeq = lpo_from_aln_mat (n_seqs, seq_names, seq_titles, aln_mat, aln_lengths, do_switch_case, score_matrix);
  
 free_memory_and_exit:
  for (i=0; i<n_seqs; i++) {
    FREE (seq_names[i]);
    FREE (seq_titles[i]);
    FREE (aln_mat[i]);
  }
  FREE (seq_names);
  FREE (seq_titles);
  FREE (aln_mat);
  FREE (aln_lengths);
  
  if (n_seqs>0) {
    strcpy (lposeq->name, lposeq->source_seq[0].name);
    FREE (lposeq->title);
    lposeq->title = strdup (lposeq->source_seq[0].title);
  }
  
  return lposeq;
}

/** Creates an LPO from an RC-MSA alignment matrix.
 */
LPOSequence_T *lpo_from_aln_mat (int n_seqs, char **seq_names, char **seq_titles, char **aln_mat, int *aln_lengths, int do_switch_case, ResidueScoreMatrix_T *score_matrix)
{
  int i, j, k, len, col, res_id, letter_id;
  char ch;
  LPOSequence_T *curr_seq, *lposeq;
  int **column_ids; /** which column contains each residue */
  int **res_ids; /** which residue is in each column */
  int *al_x, *al_y;
  int max_aln_length = 0;
  char *consens_row;

  if (n_seqs==0)
    return NULL;
  
  CALLOC (column_ids, n_seqs+1, int *);
  CALLOC (res_ids, n_seqs+1, int *);
  for (i=0; i<n_seqs; i++) {
    if (aln_lengths[i] > max_aln_length) {
      max_aln_length = aln_lengths[i];
    }
  }
  for (i=0; i<n_seqs+1; i++) {
    CALLOC (column_ids[i], max_aln_length, int);
    CALLOC (res_ids[i], max_aln_length, int);
  }

  /* GET CONSENSUS ROW THAT ALIGNS TO SOME OTHER ROW IN EACH COLUMN */
  /* CONSENSUS ROW IS _NECESSARY_ TO DEAL WITH THE POSSIBILITY THAT, e.g., UNORDERED RESIDUES
     IN seqs 0 & 1 ARE ORDERED BY ALIGNMENT WITH seq 2. */
  
  CALLOC (consens_row, max_aln_length, char);
  for (col=0; col<max_aln_length; col++) {
    consens_row[col] = '-';
    for (i=0; i<n_seqs; i++) {
      if (col < aln_lengths[i] && is_residue_char(aln_mat[i][col])) {
	consens_row[col] = aln_mat[i][col];
	break;
      }
    }
  }

  /* MAKE SEQUENCE FROM CONSENSUS ROW */

  CALLOC (lposeq, 1, LPOSequence_T);
  strcpy (lposeq->name, "consens_row");
  lposeq->title = strdup ("");
  
  CALLOC (lposeq->sequence, max_aln_length+1, char);
  for (len=col=0; col<max_aln_length; col++) {
    ch = consens_row[col];
    if (do_switch_case == switch_case_to_lower) {
      ch = tolower(ch);
    }
    else if (do_switch_case == switch_case_to_upper) {
      ch = toupper(ch);
    }
    if (is_residue_char(ch)) {
      lposeq->sequence[len] = ch;
      column_ids[0][len] = col;
      res_ids[0][col] = len;
      len++;
    }
  }
  lposeq->length = len;
  lposeq->sequence[col] = '\0';
  REALLOC (lposeq->sequence, len+1, char);
    
  initialize_seqs_as_lpo (1, lposeq, score_matrix);
  
  for (i=0; i<n_seqs; i++) {

    CALLOC (curr_seq, 1, LPOSequence_T);
    strcpy (curr_seq->name, seq_names[i]);
    curr_seq->title = strdup (seq_titles[i]);

    /* READ CHARACTERS FROM ALIGNMENT ROW INTO SEQUENCE: */
    CALLOC (curr_seq->sequence, aln_lengths[i]+1, char);
    for (len=col=0; col<aln_lengths[i]; col++) {
      ch = aln_mat[i][col];
      if (do_switch_case == switch_case_to_lower) {
	ch = tolower(ch);
      }
      else if (do_switch_case == switch_case_to_upper) {
	ch = toupper(ch);
      }
      if (is_residue_char(ch)) {
	curr_seq->sequence[len] = ch;
	column_ids[i+1][len] = col;
	res_ids[i+1][col] = len;
	len++;
      }
    }
    curr_seq->length = len;
    curr_seq->sequence[len] = '\0';
    REALLOC (curr_seq->sequence, len+1, char);
    
    initialize_seqs_as_lpo (1, curr_seq, score_matrix);
    
    /* ALIGN THIS SEQUENCE TO EXISTING ALIGNMENT USING CONSENSUS ROW */

    build_seq_to_po_index (lposeq);
    
    CALLOC (al_x, len, int);
    CALLOC (al_y, lposeq->length, int);
    
    for (j=0; j<lposeq->length; j++) {
      al_y[j] = INVALID_LETTER_POSITION;
    }
    
    /* FOR EACH RESIDUE IN THIS SEQUENCE: */
    for (j=0; j<len; j++) {
      col = column_ids[i+1][j];
      
      /* GET AN OLD NODE IN CONSENSUS ROW ALIGNED TO THIS ONE: */
      res_id = res_ids[0][col];
      
      /* STORE INFO IN al_x,al_y: */
      letter_id = lposeq->source_seq[0].seq_to_po[res_id];
      al_x[j] = letter_id;
      al_y[letter_id] = j;
    }
        
    /* DO THE FUSION-BASED BUILDUP */
    
    fuse_ring_identities (lposeq->length, lposeq->letter,
			  curr_seq->length, curr_seq->letter,
			  al_y, al_x);
    
    fuse_lpo (lposeq, curr_seq, al_y, al_x);
    
    free_lpo_sequence (curr_seq, 0);
    FREE (al_x);
    FREE (al_y);
    FREE (curr_seq);
  }
  
  /* STRIP OUT CONSENSUS SEQUENCE */
  
  strip_seq_from_lpo (lposeq, /*remove #*/ 0, /*0=don't keep all links*/ 0);
  
  for (i=0; i<n_seqs+1; i++) {
    FREE (column_ids[i]);
    FREE (res_ids[i]);
  }
  FREE (column_ids);
  FREE (res_ids);
  FREE (consens_row);
  
  return lposeq;
}


static int is_residue_char (char ch)
{
  if (ch>='a' && ch<='z') return 1;
  if (ch>='A' && ch<='Z') return 1;
  if (ch=='?' || ch=='[' || ch==']') return 1;
  return 0;
}

static int is_gap_char (char ch)
{
  if (ch=='.' || ch=='-') return 1;
  return 0;
}

static int is_name_first_char (char ch)
{
  if (ch=='\n' || ch=='\r' || ch==' ' || ch=='\t' || ch=='#' || ch=='*') return 0;
  return 1;
}

static int filter_sequence_set (int n_seqs, FILE *select_ifile, int remove_listed_sequences,
				char **seq_names, char **seq_titles, char **aln_mat, int *aln_lengths)
{
  int i, j;
  int *keep_seq;
  char name[512]="";
  
  /* DECIDE WHICH SEQUENCES TO KEEP: */
  
  CALLOC (keep_seq, n_seqs, int);
  
  for (i=0; i<n_seqs; i++) {
    keep_seq[i] = ((!select_ifile || remove_listed_sequences) ? 1 : 0);
  }
  
  while (select_ifile && fscanf (select_ifile, "SOURCENAME=%[^\r\n]%*[\r\n]", name) >= 1) {
    for (i=0; i<n_seqs; i++) {
      if (0 == strcmp(seq_names[i],name)) {
	keep_seq[i] = (remove_listed_sequences ? 0 : 1);
      }
    }
  }
  
  /* FREE MEM FOR DISCARDED SEQUENCES: */
  for (i=0; i<n_seqs; i++) if (keep_seq[i]==0) {
    FREE (seq_names[i]);
    FREE (seq_titles[i]);
    FREE (aln_mat[i]);
  }

  /* COMPACT SEQUENCE LIST (MOVE i<--j): */
  for (i=j=0; j<n_seqs; i++, j++) {
    while (j<n_seqs && keep_seq[j]==0)
      j++;
    if (j==n_seqs)
      break;
    if (j<n_seqs && i!=j) {
      seq_names[i] = seq_names[j];
      seq_titles[i] = seq_titles[j];
      aln_mat[i] = aln_mat[j];
      aln_lengths[i] = aln_lengths[j];
    }
  }
  
  FREE (keep_seq);
  
  return i;
}



void strip_seq_from_lpo (LPOSequence_T *lposeq, int bad_seq_id, int keep_all_links)
{
  int i, j, id_left, id_right, len, n_seqs;
  LPOLetterLink_T *lnk, *tmp_lnk;
  LPOLetterSource_T *src, *prev, *tmp_src;
  LPOLetter_T *lett;
  
  n_seqs = lposeq->nsource_seq - 1;
  
  /* FREE MEM ASSOCIATED WITH BAD SEQ */
  FREE (lposeq->source_seq[bad_seq_id].title);
  FREE (lposeq->source_seq[bad_seq_id].sequence);
  FREE (lposeq->source_seq[bad_seq_id].seq_to_po);
  FREE (lposeq->source_seq[bad_seq_id].po_to_seq);

  /* COMPACT SOURCE LIST */
  for (i=0; i<n_seqs; i++) if (i>=bad_seq_id) {
    lposeq->source_seq[i] = lposeq->source_seq[i+1];
  }
  /* This de-allocation causes a bug in conjunction with
     later REBUFF calls: */
  /* REALLOC (lposeq->source_seq, n_seqs, LPOSourceInfo_T); */
  lposeq->nsource_seq = n_seqs;
  
  
  /* RENUMBER SOURCES IN EACH LETTER */
  for (i=0; i<lposeq->length; i++) {
    lett = &(lposeq->letter[i]);
    prev = NULL;
    src = &(lett->source);
    while (src != NULL && src->iseq >= 0) {
      if (src->iseq == bad_seq_id) {
	if (prev) {  /* NOT IN LIST HEAD, SO RELINK FROM PREV AND FREE */
	  prev->more = src->more;
	  FREE (src);
	  src = prev->more;
	}
	else {  /* IN LIST HEAD, SO REASSIGN HEAD DATA, RELINK, FREE */
	  if (src->more) {
	    tmp_src = src->more;
	    src->ipos = tmp_src->ipos;
	    src->iseq = tmp_src->iseq;
	    src->more = tmp_src->more;
	    FREE (tmp_src);
	  }
	  else {
	    src->ipos = -1;
	    src->iseq = -1;
	  }
	}
      }
      else {
	if (src->iseq>=bad_seq_id)
	  src->iseq--;
	prev = src;
	src = src->more;
      }
    }
  }
  
  if (0==keep_all_links) {

    build_seq_to_po_index (lposeq);
    
    /* FREE ALL LINK INFO */
    for (i=0; i<lposeq->length; i++) {
      lett = &(lposeq->letter[i]);
      for (lnk = &(lett->left); lnk->more != NULL; ) {
	tmp_lnk = lnk->more;
	lnk->more = tmp_lnk->more;
	FREE (tmp_lnk);
      }
      lnk->ipos = -1;
      for (lnk = &(lett->right); lnk->more != NULL; ) {
	tmp_lnk = lnk->more;
	lnk->more = tmp_lnk->more;
	FREE (tmp_lnk);
      }
      lnk->ipos = -1;
    }
    
    /* REBUILD LINK INFO FROM REMAINING SEQS */
    for (i=0; i<n_seqs; i++) {
      len = lposeq->source_seq[i].length;
      for (j=0; j<len-1; j++) {
	id_left = lposeq->source_seq[i].seq_to_po[j];
	id_right = lposeq->source_seq[i].seq_to_po[j+1];
	add_lpo_link (&(lposeq->letter[id_left].right), id_right);
	add_lpo_link (&(lposeq->letter[id_right].left), id_left);
      }
    }
  }
}

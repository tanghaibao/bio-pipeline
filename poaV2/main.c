

#include "lpo.h"
#include "msa_format.h"
#include "align_score.h"


static LPOSequence_T *read_partial_order_file (char *po_filename, char *subset_filename, int remove_listed_seqs, int keep_all_links, int do_switch_case, ResidueScoreMatrix_T *mat);

int main(int argc,char *argv[])
{
  int i,j,ibundle=ALL_BUNDLES,nframe_seq=0,use_reverse_complement=0;
  int nseq=0,do_switch_case=dont_switch_case,do_analyze_bundles=0;
  int nseq_in_list=0,n_input_seqs=0,max_input_seqs=0;
  char score_file[256],seq_file[256],po_list_entry_filename[256],*comment=NULL,*al_name="test align";
  ResidueScoreMatrix_T score_matrix; /* DEFAULT GAP PENALTIES*/
  LPOSequence_T *seq=NULL,*lpo_out=NULL,*frame_seq=NULL,*dna_lpo=NULL,*lpo_in=NULL;
  LPOSequence_T **input_seqs=NULL;
  FILE *errfile=stderr,*logfile=NULL,*lpo_file_out=NULL,*po_list_file=NULL,*seq_ifile=NULL;
  char *print_matrix_letters=NULL,*fasta_out=NULL,*po_out=NULL,*matrix_filename=NULL,
    *seq_filename=NULL,*frame_dna_filename=NULL,*po_filename=NULL,*po2_filename=NULL,
    *po_list_filename=NULL, *hbmin=NULL,*numeric_data=NULL,*numeric_data_name="Nmiscall",
    *dna_to_aa=NULL,*pair_score_file=NULL,*aafreq_file=NULL,*termval_file=NULL,
    *bold_seq_name=NULL,*subset_file=NULL,*subset2_file=NULL,*rm_subset_file=NULL,
    *rm_subset2_file=NULL;
  float bundling_threshold=0.9;
  int exit_code=0,count_sequence_errors=0,please_print_snps=0,
    report_consensus_seqs=0,report_major_allele=0,use_aggressive_fusion=0;
  int show_allele_evidence=0,please_collapse_lines=0,keep_all_links=0;
  int remove_listed_seqs=0,remove_listed_seqs2=0,please_report_similarity;
  int do_global=0, do_progressive=0, do_preserve_sequence_order=0;
  char *reference_seq_name="CONSENS%d",*clustal_out=NULL;

  black_flag_init(argv[0],PROGRAM_VERSION);

  if (argc<2) {
    fprintf(stderr,"\nUsage: %s [OPTIONS] MATRIXFILE\n"
"Align a set of sequences or alignments using the scores in MATRIXFILE.\n"
"Example: %s -read_fasta multidom.seq -clustal m.aln blosum80.mat\n\n"
"INPUT:\n"
"  -read_fasta FILE       Read in FASTA sequence file.\n"
"  -read_msa FILE         Read in MSA alignment file.\n"
"  -read_msa2 FILE        Read in second MSA file. \n" 
"  -subset FILE           Filter MSA to include list of seqs in file.\n"
"  -subset2 FILE          Filter second MSA to include list of seqs in file.\n"
"  -remove FILE           Filter MSA to exclude list of seqs in file.\n"
"  -remove2 FILE          Filter second MSA to exclude list of seqs in file.\n"
"  -read_msa_list FILE    Read an MSA from each filename listed in file.\n"
"  -tolower               Force FASTA/MSA sequences to lowercase\n"
"                           (nucleotides in our matrix files)\n"
"  -toupper               Force FASTA/MSA sequences to UPPERCASE\n"
"                           (amino acids in our matrix files)\n"
"\nALIGNMENT:\n"
"  -do_global             Do global alignment.\n"
"  -do_progressive        Perform progressive alignment using a guide tree\n"
"                           built by neighbor joining from a set of\n"
"                           sequence-sequence similarity scores.\n"
"  -read_pairscores FILE  Read tab-delimited file of similarity scores.\n"
"                           (If not provided, scores are constructed\n"
"                           using pairwise sequence alignment.)\n"
"  -fuse_all              Fuse identical letters on align rings.\n"
"\nANALYSIS:\n"
"  -hb                    Perform heaviest bundling to generate consensi.\n"
"  -hbmin VALUE           Include in heaviest bundle sequences with\n"
"                           percent ID (as a fraction) >= value.\n"
"\nOUTPUT:\n"
"  -pir FILE              Write out MSA in PIR format.\n"
"  -clustal FILE          Write out MSA in CLUSTAL format.\n"
"  -po FILE               Write out MSA in PO format.\n"
"  -preserve_seqorder     Write out MSA with sequences in their input order.\n"
"  -printmatrix LETTERS   Print score matrix to stdout.\n"
"  -best                  Restrict MSA output to heaviest bundles (PIR only).\n"
"  -v                     Run in verbose mode (e.g. output gap penalties).\n\n"
"  NOTE:  One of the -read_fasta, -read_msa, or -read_msa_list arguments\n" 
"         must be used, since a sequence or alignment file is required.\n\n"
"For more information, see http://www.bioinformatics.ucla.edu/poa.\n\n"
	    ,argv[0],argv[0]);
    exit(-1);
  }

  FOR_ARGS(i,argc) { /* READ ALL THE ARGUMENTS */
    ARGMATCH_VAL("-tolower",do_switch_case,switch_case_to_lower);
    ARGMATCH_VAL("-toupper",do_switch_case,switch_case_to_upper);
    ARGMATCH_VAL("-v",logfile,stdout);
    ARGMATCH_VAL("-best",ibundle,0); /*RESTRICT FASTA OUTPUT TO HB */
    ARGMATCH_VAL("-hb",do_analyze_bundles,1);/*CALCULATE HEAVIEST BUNDLING*/
    ARGGET("-printmatrix",print_matrix_letters);
    ARGGET("-read_msa",po_filename); /* READ A MSA FILE FOR ALIGNMENT/ANALYSIS*/
    ARGGET("-read_msa2",po2_filename); /* READ A SECOND MSA FILE FOR ALIGNMENT/ANALYSIS*/  
    ARGGET("-read_msa_list",po_list_filename); /* READ A LIST OF MSAs FOR ALIGNMENT/ANALYSIS */
    ARGGET("-pir",fasta_out); /* SAVE FASTA-PIR FORMAT ALIGNMENT FILE */
    ARGGET("-clustal",clustal_out); /* SAVE CLUSTAL FORMAT ALIGNMENT FILE */
    ARGGET("-po",po_out); /* SAVE PO FORMAT ALIGNMENT FILE */
    ARGMATCH("-preserve_seqorder",do_preserve_sequence_order);  /* DO PRESERVE SEQUENCE ORDER */
    ARGGET("-hbmin",hbmin); /* SET THRESHOLD FOR BUNDLING */
    ARGMATCH("-fuse_all",use_aggressive_fusion);
    ARGMATCH("-do_global",do_global); /* DO GLOBAL */
    ARGGET("-read_pairscores",pair_score_file); /* FILENAME TO READ PAIR SCORES*/
    ARGMATCH("-do_progressive", do_progressive); /* DO PROGRESSIVE ALIGNMENT */
    ARGGET("-subset",subset_file); /* FILENAME TO READ SEQ SUBSET LIST*/
    ARGGET("-subset2",subset2_file); /* FILENAME TO READ SEQ SUBSET LIST*/
    ARGGET("-remove",rm_subset_file); /* FILENAME TO READ SEQ REMOVAL LIST*/
    ARGGET("-remove2",rm_subset2_file); /* FILENAME TO READ SEQ REMOVAL LIST*/
    ARGGET("-read_fasta",seq_filename); /* READ FASTA FILE FOR ALIGNMENT */
    NEXTARG(matrix_filename); /* NON-FLAG ARG SHOULD BE MATRIX FILE */   
  }

  
  /** CHECK FOR CONFLICTING FLAGS **/
  
  if (po_list_filename && (po_filename || po2_filename)) {
    WARN_MSG(USERR,(ERRTXT, "Error: The -read_po_list and -read_po flags cannot be used at the same time.\nExiting."), "$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if (((subset_file || rm_subset_file) && !po_filename) || ((subset2_file || rm_subset2_file) && !po2_filename)) {
    WARN_MSG(USERR,(ERRTXT, "Error: Each -subset/-remove flag must have a corresponding -read_po flag.\nExiting."),"$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if ((subset_file && rm_subset_file) || (subset2_file && rm_subset2_file)) {
    WARN_MSG(USERR,(ERRTXT, "Error: The -subset and -remove flags cannot be used at the same time.\nExiting."),"$Revision: 1.2.2.9 $");
    exit_code = 1;
    goto free_memory_and_exit;
  }
  
  if (rm_subset_file) {
    subset_file = rm_subset_file;
    remove_listed_seqs = 1;
  }
  if (rm_subset2_file) {
    subset2_file = rm_subset2_file;
    remove_listed_seqs2 = 1;
  }
  
  if (hbmin)
    bundling_threshold=atof(hbmin);  

  if (!matrix_filename ||
      read_score_matrix(matrix_filename,&score_matrix)<=0){/* READ MATRIX */
    WARN_MSG(USERR,(ERRTXT,"Error reading matrix file %s.\nExiting",
		    matrix_filename ? matrix_filename: "because none specified"),"$Revision: 1.2.2.9 $");
    exit_code=1; /* SIGNAL ERROR CONDITION */
    goto free_memory_and_exit;
  }
  
  if (logfile) {
    fprintf(logfile,"X-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): %d %d %d %d %d\n",
	    score_matrix.gap_penalty_set[0][0],
	    score_matrix.gap_penalty_set[0][1],
	    score_matrix.gap_penalty_set[0][2],
	    score_matrix.trunc_gap_length,
	    score_matrix.decay_gap_length);
    fprintf(logfile,"X-Gap Penalties (0, 1, 2, ...): ");
    for (i=0; i<=score_matrix.max_gap_length; i++) {
      fprintf (logfile, "%d ", score_matrix.gap_penalty_x[i]);
    }
    fprintf(logfile,"... \n");
    fprintf(logfile,"Y-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): %d %d %d %d %d\n",
	    score_matrix.gap_penalty_set[1][0],
	    score_matrix.gap_penalty_set[1][1],
	    score_matrix.gap_penalty_set[1][2],
	    score_matrix.trunc_gap_length,
	    score_matrix.decay_gap_length);
    fprintf(logfile,"Y-Gap Penalties (0, 1, 2, ...): ");
    for (i=0; i<=score_matrix.max_gap_length; i++) {
      fprintf (logfile, "%d ", score_matrix.gap_penalty_y[i]);
    }
    fprintf(logfile,"... \n");
  }
  
  if (print_matrix_letters) /* USER WANTS US TO PRINT A MATRIX */
    print_score_matrix(stdout,&score_matrix,print_matrix_letters
		       /*"ARNDCQEGHILKMFPSTWYV"*/);
  
  
  /** READ INPUT FILES **/
  
  n_input_seqs = 0;
  max_input_seqs = 10;
  CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
  
  if (po_filename) {
    lpo_in = read_partial_order_file (po_filename, subset_file, remove_listed_seqs, keep_all_links, do_switch_case, &score_matrix);
    if (lpo_in == NULL) {
      exit_code = 1;
      goto free_memory_and_exit;
    }
    fprintf(errfile,"...Read %d sequences from MSA file %s...\n",lpo_in->nsource_seq,po_filename);
    input_seqs[n_input_seqs++] = lpo_in;
    lpo_in = NULL;
  }
  
  if (po2_filename) {
    lpo_in = read_partial_order_file (po2_filename, subset2_file, remove_listed_seqs2, keep_all_links, do_switch_case, &score_matrix);
    if (lpo_in == NULL) {
      exit_code = 1;
      goto free_memory_and_exit;
    }
    fprintf(errfile,"...Read %d sequences from second MSA file %s...\n",lpo_in->nsource_seq,po2_filename);
    input_seqs[n_input_seqs++] = lpo_in;
    lpo_in = NULL;
  }
  
  if (po_list_filename) {
    po_list_file = fopen (po_list_filename, "r");
    while (po_list_file && fscanf (po_list_file, " %s", po_list_entry_filename) == 1) {
      lpo_in = read_partial_order_file (po_list_entry_filename, NULL, 0, 0, do_switch_case, &score_matrix);
      if (lpo_in == NULL) {
	exit_code = 1;
	goto free_memory_and_exit;
      }
      fprintf(errfile,"...Read %d sequences from PO list entry %s...\n",lpo_in->nsource_seq,po_list_entry_filename);
      nseq_in_list += lpo_in->nsource_seq;
      input_seqs[n_input_seqs++] = lpo_in;
      lpo_in = NULL;
      if (n_input_seqs == max_input_seqs) {
	max_input_seqs *= 2;
	REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
      }
    }
    if (nseq_in_list==0) {
      WARN_MSG(USERR,(ERRTXT,"Error reading PO list file %s.\nExiting",
		      po_list_file),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */
      goto free_memory_and_exit;      
    }
  }
  
  if (seq_filename) {
    seq_ifile = fopen (seq_filename, "r");
    if (seq_ifile == NULL) {
      WARN_MSG(USERR,(ERRTXT,"Couldn't open sequence file %s.\nExiting",
		      seq_filename),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */
      goto free_memory_and_exit;
    }
    nseq = read_fasta (seq_ifile, &seq, do_switch_case, &comment);
    fclose (seq_ifile);
    if (nseq == 0) {
      WARN_MSG(USERR,(ERRTXT,"Error reading sequence file %s.\nExiting",
		      seq_filename),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */
      goto free_memory_and_exit;
    }
    fprintf(errfile,"...Read %d sequences from sequence file %s...\n",nseq,seq_filename);
    for (i=0; i<nseq; i++) {
      input_seqs[n_input_seqs++] = &(seq[i]);
      initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
      if (n_input_seqs == max_input_seqs) {
	max_input_seqs *= 2;
	REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
      }
    }
  }


  /** BUILD AND ANALYZE OUTPUT PO-MSA **/

  if (n_input_seqs == 0) { /* HMM.. NO DATA. */
    WARN_MSG(USERR,(ERRTXT,"No input sequences were provided; use one of the -read_ flags.\nExiting."), "$Revision: 1.2.2.9 $");
    exit_code=1; /* SIGNAL ERROR CONDITION */
    goto free_memory_and_exit;
  }
  else {
    lpo_out = buildup_progressive_lpo (n_input_seqs, input_seqs, &score_matrix,
				       use_aggressive_fusion, do_progressive, pair_score_file,
				       matrix_scoring_function, do_global, do_preserve_sequence_order);
  }
  
  if (comment) { /* SAVE THE COMMENT LINE AS TITLE OF OUR LPO */
    FREE(lpo_out->title);
    lpo_out->title=strdup(comment);
  }
  
  /* DIVIDE INTO BUNDLES W/ CONSENSUS USING PERCENT ID */
  if (do_analyze_bundles)
    generate_lpo_bundles(lpo_out,bundling_threshold);

  if (po_out) { /* WRITE FINAL PARTIAL ORDER ALIGNMENT TO OUTPUT */
    if (lpo_file_out=fopen(po_out, "w")) {
       write_lpo(lpo_file_out,lpo_out,&score_matrix);
       fclose(lpo_file_out);
       fprintf(errfile,"...Wrote %d sequences to PO file %s...\n",lpo_out->nsource_seq,po_out);
    }
    else {
      WARN_MSG(USERR,(ERRTXT,"*** Could not save PO file %s.  Exiting.",
		      po_out),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */ 
    }
  }

  if (fasta_out) { /* WRITE FINAL ALIGNMENT IN FASTA-PIR FORMAT */
    if (seq_ifile=fopen(fasta_out,"w")) { /* FASTA-PIR ALIGNMENT*/
      write_lpo_bundle_as_fasta(seq_ifile,lpo_out,score_matrix.nsymbol,
				score_matrix.symbol,ibundle);
      fclose(seq_ifile);
      fprintf(errfile,"...Wrote %d sequences to FASTA-PIR file %s...\n",lpo_out->nsource_seq,fasta_out);
    }
    else {
      WARN_MSG(USERR,(ERRTXT,"*** Could not save FASTA-PIR file %s.  Exiting.",
		      fasta_out),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */ 
   }
  }

  if (clustal_out) { /* WRITE FINAL ALIGNMENT IN CLUSTAL FORMAT */
    if (seq_ifile=fopen(clustal_out,"w")) { /* CLUSTAL ALIGNMENT*/
      export_clustal_seqal(seq_ifile,lpo_out,score_matrix.nsymbol,
			   score_matrix.symbol);
      fclose(seq_ifile);
      fprintf(errfile,"...Wrote %d sequences to CLUSTAL file %s...\n",lpo_out->nsource_seq,clustal_out);
    }
    else {
      WARN_MSG(USERR,(ERRTXT,"*** Could not save CLUSTAL file %s.  Exiting.",
		      clustal_out),"$Revision: 1.2.2.9 $");
      exit_code=1; /* SIGNAL ERROR CONDITION */ 
   }
  }



 free_memory_and_exit: /* FREE ALL DYNAMICALLY ALLOCATED DATA!!!! */
  
  if (dna_lpo)
    free_lpo_sequence(dna_lpo,TRUE);
  
  for (i=0;i<n_input_seqs;i++) {
    for (j=0;j<nseq;j++) {
      if (input_seqs[i]==&(seq[j]))
	break;
    }
    free_lpo_sequence(input_seqs[i],(j==nseq));
  }
  FREE (input_seqs);
  if (nseq>0) FREE (seq);
  
  exit (exit_code);
}



static LPOSequence_T *read_partial_order_file (char *po_filename, char *subset_filename, int remove_listed_seqs, int keep_all_links, int do_switch_case, ResidueScoreMatrix_T *mat)
{
  LPOSequence_T *lpo_in;
  FILE *po_file=NULL, *subset_file=NULL;
  
  if (!po_filename)
    return NULL;
  
  po_file = fopen (po_filename, "r");
  if (!po_file) {
    WARN_MSG (USERR, (ERRTXT,"Couldn't open MSA file %s.\nExiting.",po_filename), "$Revision: 1.2.2.9 $");
    return NULL;
  }
  
  if (subset_filename) {
    subset_file = fopen (subset_filename, "r");
    if (!subset_file) {
      WARN_MSG (USERR, (ERRTXT,"Couldn't open subset file %s.\nExiting.",subset_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }
  
  if (subset_file) {
    lpo_in = read_msa_select (po_file, UNKNOWN_MSA, subset_file, keep_all_links, remove_listed_seqs, do_switch_case, mat);
    fclose (subset_file);
    fclose (po_file);
    if (lpo_in==NULL || lpo_in->nsource_seq == 0) {
      WARN_MSG (USERR, (ERRTXT,"MSA file %s, filtered with subset file %s, couldn't be read or contains no sequences.\nExiting.", po_filename, subset_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }
  else {
    lpo_in = read_msa (po_file, UNKNOWN_MSA, do_switch_case, mat);
    fclose (po_file);
    if (lpo_in==NULL || lpo_in->nsource_seq == 0) {
      WARN_MSG (USERR, (ERRTXT,"MSA file %s couldn't be read or contains no sequences.\nExiting.", po_filename), "$Revision: 1.2.2.9 $");
      return NULL;
    }
  }

  return lpo_in;
}

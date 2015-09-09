

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/** if two align-rings are aligned to each other, make sure
    that the (single) aligned residue pair consists of identical
    residues, if possible.
    :
    ((a,-),(b,-),(c,d),(-,a),(-b)) ==>
    ((a,a),(b,-),(c,-),(-,d),(-b)) OR ((a,-),(b,b),(c,-),(-,d)).
    :
    ((a,-),(c,d),(-,b)) ==> self.
*/
void fuse_ring_identities(int len_x,LPOLetter_T seq_x[],
			  int len_y,LPOLetter_T seq_y[],
			  LPOLetterRef_T al_x[],
			  LPOLetterRef_T al_y[])
{
  int i,j;
  LOOP (i,len_y) {
    if (al_y[i]<0 || seq_x[al_y[i]].letter == seq_y[i].letter)
      continue; /* NOT ALIGNED, OR ALREADY IDENTICAL, SO SKIP */
    for (j=seq_x[al_y[i]].align_ring;j!=al_y[i];j=seq_x[j].align_ring) {
      if (seq_x[j].letter == seq_y[i].letter) { /* IDENTICAL! SO FUSE! */
	al_x[al_y[i]]= INVALID_LETTER_POSITION; /* DISCONNECT FROM OLD */
	al_y[i]=j; /* CONNECT TO NEW IDENTITY */
	al_x[j]=i;
	break; /* SEARCH YE NO FURTHER */
      }
    }
  }
}

/** if two align-rings are aligned to each other, make sure
    that as many identical-residue pairs are fused as possible.
    :
    ((a,-),(b,-),(c,d),(-,a),(-b)) ==>
    ((a,a),(b,b),(c,-),(-,d)).
    :
    ((a,-),(c,d),(-,b)) ==> self.

    NB:
    THIS DOES NOT WORK with the current fuse_lpo function.
    The fusion 
*/
void full_fuse_ring_identities (int len_x, LPOLetter_T *seq_x,
				int len_y, LPOLetter_T *seq_y,
				LPOLetterRef_T *al_x,
				LPOLetterRef_T *al_y)
{
  int i, ip, j, jp;  /* i LABELS POS IN seq_x, j LABELS POS IN seq_y */
  int ck;  /* WAS IDENTICAL-RESIDUE PAIR FOUND? */
  for (i=0; i<len_x; i++) if ((j = al_x[i]) >= 0) {
    al_x[i] = al_y[j] = INVALID_LETTER_POSITION;  /* DISCONNECT FROM OLD */
    
    /* i,j ARE AN ALIGNED PAIR.  WALK THROUGH RESPECTIVE RINGS: */
    ip=i; jp=j; ck=0;
    do /* while (ip!=i) */ {
      do /* while (jp!=j) */ {
	if (seq_x[ip].letter == seq_y[jp].letter) {  /* IDENTICAL! SO FUSE! */
	  al_x[ip] = jp;
	  al_y[jp] = ip;
	  ck=1;
	  jp=j;  /* EXIT TO OUTER LOOP... AT MOST ONE FUSED TO EACH POS IN seq_y. */
	}
	else {
	  jp = seq_y[jp].align_ring;
	}
      }
      while (jp!=j);
      ip = seq_x[ip].align_ring;
    }
    while (ip!=i);
    
    if (ck==0) {  /* NO IDENTICAL-RESIDUE PAIR FOUND, SO RECONNECT ORIGINAL */
      al_x[i] = j;
      al_y[j] = i;
    }
  }
}


/** aligns the sequences in seq[] to the sequence or partial order in
  new_seq; seq[] must be linear orders (regular sequences);
  the alignment is built up by iterative partial order alignment,
  and the resulting partial order is returned in new_seq */
LPOSequence_T *buildup_lpo(LPOSequence_T *new_seq,
			   int nseq,LPOSequence_T seq[],
			   ResidueScoreMatrix_T *score_matrix,
			   int use_aggressive_fusion,
                           int use_global_alignment)
{
  int i,max_alloc=0,total_alloc;
  LPOLetterRef_T *al1=NULL,*al2=NULL;

  lpo_index_symbols(new_seq,score_matrix); /* MAKE SURE LPO IS TRANSLATED */
  for (i=0;i<nseq;i++) { /* ALIGN ALL SEQUENCES TO my_lpo ONE BY ONE */
    if (seq[i].letter == NULL) /* HMM.  HASN'T BEEN INITIALIZED AT ALL YET */
      initialize_seqs_as_lpo(1,seq+i,score_matrix);
    total_alloc=new_seq->length*seq[i].length
      + sizeof(LPOLetter_T)*new_seq->length;
    if (total_alloc>max_alloc) { /* DP RECTANGLE ARRAY SIZE */
      max_alloc=total_alloc;
#ifdef REPORT_MAX_ALLOC
      fprintf(stderr,"max_alloc: %d bytes\n",max_alloc);
#endif
      if (max_alloc>POA_MAX_ALLOC) {
	WARN_MSG(TRAP,(ERRTXT,"Exceeded memory bound: %d\n Exiting!\n\n",max_alloc),"$Revision: 1.2.2.9 $");
	break; /* JUST RETURN AND FINISH */
      }
    }
    align_lpo_po (new_seq,&seq[i],
		  score_matrix,&al1,&al2,NULL,use_global_alignment); /* ALIGN ONE MORE SEQ */
    if (use_aggressive_fusion) 
      fuse_ring_identities(new_seq->length,new_seq->letter,
			   seq[i].length,seq[i].letter,al1,al2);
    fuse_lpo(new_seq,seq+i,al1,al2); /* BUILD COMPOSITE LPO */

    free_lpo_letters(seq[i].length,seq[i].letter,TRUE);/*NO NEED TO KEEP*/
    seq[i].letter=NULL; /* MARK AS FREED... DON'T LEAVE DANGLING POINTER! */
    FREE(al1); /* DUMP TEMPORARY MAPPING ARRAYS */
    FREE(al2);
  }

  return new_seq;
}
/**@memo example: aligning a set of sequences to a partial order: 
      lpo_out=buildup_lpo(lpo_in,nseq,seq,&score_matrix);
*/



/** CLIPS seq->letter[] TO JUST THE SEGMENT ALIGNED TO letter_x[] via al_x[]
 DOES *NOT* FREE existing seq->letter[]; YOU MUST KEEP IT OR FREE IT YOURSELF*/
int clip_unaligned_ends(LPOSequence_T *seq,
			LPOLetterRef_T al[],
			int len_x,LPOLetter_T letter_x[],
			LPOLetterRef_T al_x[],int *offset,int *match_length)
{
  int i,j=0,start,end,new_length,allow_end_length=0,nidentity=0;
  LPOLetter_T *temp=NULL;
  CALLOC(temp,seq->length,LPOLetter_T); /* ALLOCATE NEW letter[] COPY */
  for (start=0;start<seq->length;start++) /* FIND 1ST ALIGNED POS */
    if (al[start]>=0)
      break;

  for (end=seq->length -1;end>=0;end--) /* FIND LAST ALIGNED POS */
    if (al[end]>=0)
      break;

  for (i=start;i<=end;i++) /* COUNT IDENTITIES TO letter_x[] */
    if (al[i]>=0 && seq->letter[i].letter==letter_x[al[i]].letter)
      nidentity++;
  if (match_length) /* RETURN THE MATCH LENGTH TO THE CALLER */
    *match_length = end-start+1;

  if (start>allow_end_length) /* ALLOW EXTRA RESIDUES ON EITHER END*/
    start-=allow_end_length;
  else /* KEEP IN BOUNDS */
    start=0;
  if (end+allow_end_length<seq->length)
    end+=allow_end_length;
  else /* KEEP IN BOUNDS */
    end=seq->length-1;

  LOOP (i,len_x) /* WE ARE SHIFTING al TO THE RIGHT BY start POSITIONS */
    if (al_x[i]>=0) /* SO WE HAVE TO TRANSLATE al_x CORRESPONDINGLY */
      al_x[i]-= start;

  seq->length=end-start+1; /* NOW TRANSLATE left, right, align_ring, ring_id*/
  memcpy(temp,seq->letter+start,sizeof(LPOLetter_T)*(seq->length));
  LOOP (i,seq->length) { /* THIS *ONLY* WORKS FOR PURE LINEAR SEQUENCE!!! */
    temp[i].left.ipos -= start; /*IF <0, BECOMES INVALID BY DEFINITION, OK*/
    temp[i].right.ipos -= start;
    if (temp[i].right.ipos>=seq->length) /* PAST THE NEW, CLIPPED END */
      temp[i].right.ipos= INVALID_LETTER_POSITION;
    temp[i].ring_id=temp[i].align_ring=i;
  }

  if (offset) /* RETURN THE OFFSET TO THE CALLER */
    *offset = start;

  seq->letter=temp; /* NEW START: FIRST ALIGNED POSITION */
  return nidentity; /* NEW LENGTH: FROM 1ST TO LAST ALIGNED POS*/
}





void restore_lpo_size(LPOSequence_T *seq,int length,LPOLetter_T *letter)
{

  free_lpo_letters(seq->length,seq->letter,TRUE); /* DUMP CLIPPED VERSION*/
  seq->length=length; /* RESTORE ORIGINAL length AND letter[] */
  seq->letter=letter;
}



/** BUILDS UP ALIGNMENT, BUT CLIPS UNALIGNED ENDS OF EACH NEW SEQUENCE ADDED  
-------------------------------------------------------
---------------------------------------------------------------------------
*/
LPOSequence_T *buildup_clipped_lpo(LPOSequence_T *new_seq,
				   int nseq,LPOSequence_T seq[],
				   ResidueScoreMatrix_T *score_matrix,
                                   int use_global_alignment)
{
  int i,ntemp,offset=0,nidentity,length_max=0,match_length=0;
  int total_alloc,max_alloc=0;
  LPOLetterRef_T *al1=NULL,*al2=NULL;
  LPOLetter_T *temp;
  float identity_max=0.,f;

  lpo_index_symbols(new_seq,score_matrix); /* MAKE SURE LPO IS TRANSLATED */
  for (i=0;i<nseq;i++) { /* ALIGN ALL SEQUENCES TO new_seq ONE BY ONE */
    if (seq[i].letter == NULL) /* HMM.  HASN'T BEEN INITIALIZED AT ALL YET */
      initialize_seqs_as_lpo(1,seq+i,score_matrix);
    total_alloc=new_seq->length*seq[i].length
      + sizeof(LPOLetter_T)*new_seq->length;
    if (total_alloc>max_alloc) { /* DP RECTANGLE ARRAY SIZE */
      max_alloc=total_alloc;
#ifdef REPORT_MAX_ALLOC
      fprintf(stderr,"max_alloc: %d bytes (%d x %d)\n",max_alloc,
	      new_seq->length,seq[i].length);
#endif
      if (max_alloc>POA_MAX_ALLOC) {
	WARN_MSG(TRAP,(ERRTXT,"Exceeded memory bound: %d\n Exiting!\n\n",max_alloc),"$Revision: 1.2.2.9 $");
	break; /* JUST RETURN AND FINISH */
      }
    }
    align_lpo_po (new_seq, &seq[i],
		  score_matrix,&al1,&al2,NULL,use_global_alignment); /* ALIGN ONE MORE SEQ */
    ntemp=seq[i].length; /* SAVE letter[] BEFORE CLIPPING IT TO ALIGNED AREA*/
    temp=seq[i].letter;
    if ((nidentity=clip_unaligned_ends(seq+i,al2,/*THERE IS AN ALIGNED REGION*/
			new_seq->length,new_seq->letter,al1,&offset,
				       &match_length))>0) {
      f=nidentity/(float)match_length; /* CALCULATE IDENTITY FRACTION */
      if (0==i /*f>identity_max*/) { /* REPORT IDENTITY OF TOP HIT */
	identity_max=nidentity;
	length_max=match_length;
      }
      fuse_lpo(new_seq,seq+i,al1,al2+offset); /*ADD CLIPPED REGION TO LPO*/
    }
    restore_lpo_size(seq+i,ntemp,temp); /* REVERT FROM CLIPPED TO ORIGINAL*/
    FREE(al1); /* DUMP TEMPORARY MAPPING ARRAYS FROM align_lpo() */
    FREE(al2);
  }

  fprintf(stderr,"%s\tmaximum identity\t%3.1f%%\t%.0f/%d\n",new_seq->name,
	  100*identity_max/length_max,identity_max,length_max);
  return new_seq;
}


/** which LPOSeq is called, or holds a sequence called, `name'? */
int find_seq_name (int nseq, LPOSequence_T **seq, char name[])
{
  int i,j;
  for (i=0;i<nseq;i++) if (seq[i]) {
    if (0==strcmp(seq[i]->name,name))
      return i;
    for (j=0;j<seq[i]->nsource_seq;j++) {
      if (0==strcmp(seq[i]->source_seq[j].name,name))
	return i;
    }
  }
  return -1;
}


typedef struct {
  double score;
  int i;
  int j;
}
SeqPairScore_T;


/* SORT IN DESCENDING ORDER BY score (SO HIGH SIMILARITY SCORES MERGE FIRST). */
/* FOR TIES, USE ITERATIVE MERGE ORDER (1-2, then 1-3, then 1-4, etc.) */
int seqpair_score_qsort_cmp (const void *void_a, const void *void_b)
{
  const SeqPairScore_T *a = (const SeqPairScore_T *)void_a;
  const SeqPairScore_T *b = (const SeqPairScore_T *)void_b;

  if (a->score > b->score)
    return -1;
  else if (a->score < b->score)
    return 1;
  
  if (a->i > b->i)
    return 1;
  else if (a->i < b->i)
    return -1;
  
  if (a->j > b->j)
    return 1;
  else if (a->j < b->j)
    return -1;
  
  return 0;
}


SeqPairScore_T *read_seqpair_scorefile (int nseq, LPOSequence_T **seq,
					ResidueScoreMatrix_T *score_matrix,
					LPOScore_T (*scoring_function)
					(int,int,LPOLetter_T [],LPOLetter_T [],
					 ResidueScoreMatrix_T *),
					int use_global_alignment,
					int do_progressive, FILE *ifile, int *p_nscore)
{
  int i,j,nscore=0,max_nscore=0;
  int *adj_score = NULL;
  SeqPairScore_T *score_list=NULL;
  LPOLetterRef_T *al1=NULL,*al2=NULL;
  double x, min_score=0.0;
  char name1[256],name2[256];
  
  CALLOC (adj_score, nseq, int);
  
  max_nscore = nseq*nseq;
  CALLOC (score_list, max_nscore, SeqPairScore_T);
    
  if (ifile) { /* IF PAIR SCORE FILE (PROGRESSIVE ASSUMED) */
    while (fscanf(ifile," %s %s %lf",name1,name2,&x)==3) {  /* READ SCORE FILE */
      i=find_seq_name(nseq,seq,name1);
      j=find_seq_name(nseq,seq,name2);
      if (i<0 || j<0) {
	WARN_MSG(USERR,(ERRTXT,"invalid sequence pair, not found: %s,%s",name1,name2),"$Revision: 1.2.2.9 $");
	FREE (score_list);
	FREE (adj_score);
	return NULL;
      }
      
      /* fprintf(stderr,"i=%d,j=%d,x=%.2f\n",i,j,x); */
      fprintf(stderr,"Saving score from file %d (%s), %d (%s) : %.2f\n",i,name1,j,name2,x);
      if (i<j) { int swap=i;i=j;j=swap; }  /* DON'T SAVE UPPER (DUPLICATE?) HALF OF THE MATRIX */
      score_list[nscore].i = i;
      score_list[nscore].j = j;
      score_list[nscore].score = x;
      if (x<min_score) min_score = x;
      nscore++;
      if (nscore==max_nscore) {
	max_nscore *= 2;
	REALLOC (score_list, max_nscore, SeqPairScore_T);
      }
      if (j==i-1) {
	adj_score[i]=1;
      }
    }
  }
  else if (do_progressive) { /* IF PROGRESSIVE BUT NO PAIR SCORE FILE */
    for (i=0;i<nseq;i++) for (j=0;j<i;j++) { /* SCORE IS BASED ON LOCAL ALIGNMENT */
      x = align_lpo_po (seq[i],seq[j],
			score_matrix,&al1,&al2,scoring_function,use_global_alignment); 
      FREE(al1); /* DUMP TEMPORARY MAPPING ARRAYS */
      FREE(al2);
      fprintf(stderr,"Saving alignment score %d (%s), %d (%s) : %.2f\n",i,seq[i]->name,j,seq[j]->name,x);
      score_list[nscore].i = i;
      score_list[nscore].j = j;
      score_list[nscore].score = x;
      if (x<min_score) min_score = x;
      nscore++;
      if (nscore==max_nscore) {
	max_nscore *= 2;
	REALLOC (score_list, max_nscore, SeqPairScore_T);
      }
      if (j==i-1) {
	adj_score[i]=1;
      }
    }
  }
  else {  /* NOT DOING PROGRESSIVE ALIGNMENT */
    /* USE DEFAULT (=0.0) PAIRSCORES, ENSURING ITERATIVE ALIGNMENT: */
    fprintf(stderr,"Performing iterative alignment...\n");
  }
  
  for (i=1;i<nseq;i++) {
    if (0 == adj_score[i]) {
      score_list[nscore].i = i;
      score_list[nscore].j = i-1;
      score_list[nscore].score = min_score - 1.0;
      nscore++;
      if (nscore==max_nscore) {
	max_nscore *= 2;
	REALLOC (score_list, max_nscore, SeqPairScore_T);
      }
    }
  }
  FREE (adj_score);
  
  /* NOW SORT LIST IN DESCENDING ORDER AND HAND BACK TO CALLER: */
  qsort(score_list,nscore,sizeof(SeqPairScore_T),seqpair_score_qsort_cmp);
  if (p_nscore) /* RETURN LENGTH OF PAIR SCORE TABLE IF REQUESTED */
    *p_nscore=nscore;
  return score_list;
}


LPOSequence_T *buildup_progressive_lpo(int nseq,LPOSequence_T **all_seqs,
				       ResidueScoreMatrix_T *score_matrix,
				       int use_aggressive_fusion,
                                       int do_progressive,
				       char score_file[], 
				       LPOScore_T (*scoring_function)
				       (int,int,LPOLetter_T [],LPOLetter_T [],
					ResidueScoreMatrix_T *),
                                       int use_global_alignment,
				       int preserve_sequence_order)
{
  int i,j,k,max_alloc=0,total_alloc,min_counts=0;
  SeqPairScore_T *score=NULL;
  LPOSequence_T *new_seq=NULL;
  FILE *ifile=NULL;
  int *seq_cluster=NULL,cluster_i,cluster_j,nscore=0,iscore;
  int *initial_nseq, *cluster_size, *seq_id_in_cluster;
  int nseq_tot;
  
  
  /* INITIALIZE ALL UNINITIALIZED SEQS: */
  for (i=0;i<nseq;i++) {
    if (all_seqs[i]->letter == NULL) {
      initialize_seqs_as_lpo(1,all_seqs[i],score_matrix);
    }
    lpo_index_symbols(all_seqs[i],score_matrix); /* MAKE SURE LPO IS TRANSLATED */  
  }

  /* RETURN IF NOTHING TO ALIGN */
  if (nseq<=0)
    return NULL;
  else if (nseq==1)
    return all_seqs[0];
  
  
  new_seq = all_seqs[0];
    
  CALLOC(seq_cluster,nseq,int);  /* MAPS SEQS (or CLUSTERS) TO CLUSTER THEY'RE IN */
  CALLOC(seq_id_in_cluster,nseq,int);  /* INDEXES SEQS (or CLUSTERS) WITHIN EACH CLUSTER */
  CALLOC(cluster_size,nseq,int);  /* COUNTS SEQS IN SAME CLUSTER (updated w/ merges) */
  CALLOC(initial_nseq,nseq,int);  /* COUNTS SEQS INITIALLY IN EACH CLUSTER (not updated w/ merges) */
  
  for (i=nseq_tot=0;i<nseq;i++) {
    seq_cluster[i] = i;  /* CREATE TRIVIAL MAPPING, EACH SEQ ITS OWN CLUSTER */
    seq_id_in_cluster[i] = 0;
    cluster_size[i] = all_seqs[i]->nsource_seq;
    initial_nseq[i] = cluster_size[i];
    nseq_tot += cluster_size[i];
  }
    
  if (score_file) {
    ifile=fopen(score_file,"r");
    if (ifile==NULL) {
      WARN_MSG(USERR,(ERRTXT,"Error reading pair score file %s.\nExiting",
			score_file),"$Revision: 1.2.2.9 $");
      goto free_and_exit;
    }
  }
  else {
    ifile = NULL;
  }
  
  
  score = read_seqpair_scorefile(nseq,all_seqs,score_matrix,scoring_function,use_global_alignment,
				 do_progressive,ifile,&nscore);
  if (score==NULL) {
    WARN_MSG(USERR,(ERRTXT,"Error generating pair scores (file %s).\nExiting",
		    score_file ? score_file : "unspecified"),"$Revision: 1.2.2.9 $");
    goto free_and_exit;
  }
  if (ifile)
    fclose (ifile);
  
  for (iscore=0;iscore<nscore;iscore++) {
    
    /* NB: NEW CLUSTER ID WILL BE MINIMUM OF INPUT IDs,
       SO MASTER CLUSTER WILL ALWAYS BE CLUSTER 0. */
    if (seq_cluster[score[iscore].i] < seq_cluster[score[iscore].j]) {
      cluster_i=seq_cluster[score[iscore].i];
      cluster_j=seq_cluster[score[iscore].j];
    }
    else if (seq_cluster[score[iscore].j] < seq_cluster[score[iscore].i]) {
      cluster_i=seq_cluster[score[iscore].j];
      cluster_j=seq_cluster[score[iscore].i];
    }
    else /* CLUSTERS ALREADY FUSED, SO SKIP THIS PAIR */
      continue;

    fprintf(stderr,"Fusing cluster %d (%s, nseq=%d) --> %d (%s, nseq=%d)... score %.2f\n",
	    cluster_j,all_seqs[cluster_j]->name,all_seqs[cluster_j]->nsource_seq,
	    cluster_i,all_seqs[cluster_i]->name,all_seqs[cluster_i]->nsource_seq,
	    score[iscore].score);
    
    new_seq = all_seqs[cluster_i];
    total_alloc = new_seq->length * (sizeof(LPOLetter_T) + all_seqs[cluster_j]->length);
    if (total_alloc>max_alloc) { /* DP RECTANGLE ARRAY SIZE */
      max_alloc=total_alloc;
#ifdef REPORT_MAX_ALLOC
      fprintf(stderr,"max_alloc: %d bytes\n",max_alloc);
#endif
      if (max_alloc>POA_MAX_ALLOC) {
	WARN_MSG(TRAP,(ERRTXT,"Exceeded memory bound: %d\n Exiting!\n\n",max_alloc),"$Revision: 1.2.2.9 $");
	break; /* JUST RETURN AND FINISH */
      }
    }
    
#ifdef USE_LOCAL_NEUTRALITY_CORRECTION /* NO LONGER USED */
    if (score_matrix->nfreq>0) { /* CALCULATE BALANCED SCORING ON EACH PO */
      balance_matrix_score(new_seq->length,new_seq->letter,score_matrix);
      balance_matrix_score(all_seqs[cluster_j]->length,all_seqs[cluster_j]->letter,
			   score_matrix);
    }
#endif

    buildup_pairwise_lpo(new_seq,all_seqs[cluster_j],score_matrix,
			 use_aggressive_fusion,
                         scoring_function,use_global_alignment);

    LOOP (i,nseq) {  /* APPEND ALL MEMBERS OF cluster_j TO cluster_i */
      if (seq_cluster[i] == cluster_j) {
	seq_cluster[i] = cluster_i;
      	seq_id_in_cluster[i] += cluster_size[cluster_i];
      }
    }
    cluster_size[cluster_i] += cluster_size[cluster_j];
    cluster_size[cluster_j] = 0;
  }
  
  if (preserve_sequence_order) {  /* PUT SEQUENCES WITHIN LPO BACK IN THEIR ORIGINAL ORDER: */
    int *perm;
    CALLOC (perm, nseq_tot, int);

    for (i=nseq_tot=0; i<nseq; i++) {
      for (j=0; j<initial_nseq[i]; j++) {
	perm[seq_id_in_cluster[i] + j] = (nseq_tot++);
      }
    }
    for (i=0; i<nseq_tot; i++) printf ("%d ", perm[i]); printf ("\n");

    reindex_lpo_source_seqs (new_seq, perm);
    FREE (perm);
  }
  
  free_and_exit:
  FREE (initial_nseq);
  FREE (seq_cluster);
  FREE (cluster_size);
  FREE (seq_id_in_cluster);
  FREE (score);
  
  return new_seq; /* RETURN THE FINAL MASTER CLUSTER... */
}


LPOSequence_T *buildup_pairwise_lpo(LPOSequence_T seq1[],LPOSequence_T seq2[],
				    ResidueScoreMatrix_T *score_matrix,
				    int use_aggressive_fusion,
       				    LPOScore_T (*scoring_function)
				    (int,int,LPOLetter_T [],LPOLetter_T [],
				     ResidueScoreMatrix_T *),
                                    int use_global_alignment)
{
  int min_counts1=0;
  int min_counts2=0;
  LPOLetterRef_T *al1=NULL,*al2=NULL;

  lpo_index_symbols(seq1,score_matrix); /* MAKE SURE LPO IS TRANSLATED */  
  lpo_index_symbols(seq2,score_matrix); /* MAKE SURE LPO IS TRANSLATED */
  align_lpo_po (seq1, seq2, score_matrix, &al1, &al2,
		scoring_function, use_global_alignment); /* ALIGN TWO POS */
  if (use_aggressive_fusion) 
     fuse_ring_identities(seq1->length,seq1->letter,
			  seq2->length,seq2->letter,al1,al2);
  fuse_lpo(seq1,seq2,al1,al2); /* BUILD COMPOSITE LPO */
  
  /* FREE LETTERS IN SECOND LPO */
  free_lpo_letters(seq2->length,seq2->letter,TRUE);
  seq2->letter=NULL; /*MARK AS FREED. DON'T LEAVE DANGLING POINTER*/
  FREE(al1); /* DUMP TEMPORARY MAPPING ARRAYS */
  FREE(al2);
  return seq1; /* RETURN THE FINAL LPO */
}


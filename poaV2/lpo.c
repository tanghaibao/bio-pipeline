

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/** INITIALIZES LINEARIZED PARTIAL ORDER DATA STRUCTURES FOR
   A LINEAR SEQUENCE */
void lpo_init(LPOSequence_T *seq)
{
  int i;

  CALLOC(seq->letter,seq->length,LPOLetter_T);

  LOOP (i,seq->length) {
    seq->letter[i].left.ipos = SEQ_Y_LEFT(i); /* JUST A LINEAR SEQ */
    seq->letter[i].right.ipos= SEQ_Y_RIGHT(i);
    seq->letter[i].source.iseq=0; /* TRIVIAL SOURCE: POINT TO SELF */
    seq->letter[i].source.ipos=i;
    seq->letter[i].align_ring = seq->letter[i].ring_id=i; /* POINT AT SELF */
    seq->letter[i].letter = seq->sequence[i]; /* COPY OUR AA LETTER */
  }
  seq->letter[seq->length -1].right.ipos= INVALID_LETTER_POSITION;
  /* NB: letter[0].left.ipos IS INVALID_LETTER_POSITION THANKS TO SEQ_Y_LEFT()
     ABOVE */
  /* BECAUSE PO CAN CONTAIN MULTIPLE SEQUENCES, WE ALSO KEEP A SOURCE LIST.
     FOR PURE LINEAR SEQUENCE, THE LIST IS JUST OUR STARTING SEQUENCE */
  save_lpo_source(seq,seq->name,seq->title,seq->length,1,NO_BUNDLE,0,NULL);
  return;
}





/**@memo initialize one or more regular sequences (linear orders) to LPO form.
  This step is REQUIRED before running partial order alignment.  
  This routine processes each sequence with limit_residues() and 
  index_symbols(), then builds a linear LPO using lpo_init(). */
  
void initialize_seqs_as_lpo(int nseq, Sequence_T seq[],ResidueScoreMatrix_T *m)
{
  int i;
  LOOP (i,nseq) {/* EXCLUDE LETTERS THAT AREN'T IN MATRIX */
    limit_residues(seq[i].sequence,m->symbol);
    /* TRANSLATE FROM ASCII LETTERS TO COMPACTED NUMBERICAL INDEX*/
    index_symbols(seq[i].length,seq[i].sequence,seq[i].sequence,
		  m->nsymbol,m->symbol);
    lpo_init(seq+i); /* CREATE TRIVIAL, LINEAR SEQUENCE LPO */
  }
}


/** translates letter symbols on lpo.letter[i] to indexes used in scoring 
    matrix*/
void lpo_index_symbols(Sequence_T *lpo,ResidueScoreMatrix_T *m)
{
  int i;

  if (lpo->letter == NULL) { /* HMM.  HASN'T BEEN INITIALIZED AT ALL YET */
    initialize_seqs_as_lpo(1,lpo,m);
    return;
  }
  if (lpo->letter[0].letter < m->nsymbol)
    return; /* LOOKS LIKE IT'S ALREADY TRANSLATED TO INDEXES */
  LOOP (i,lpo->length) /* READ FROM FILE, MAY NOT BE TRANSLATED YET */
    index_symbols(1,&lpo->letter[i].letter,&lpo->letter[i].letter,
		  m->nsymbol,m->symbol);
}





/** finds ipos for the designated sequence in the designated letter,
    or INVALID_LETTER_POSITION if iseq is not found. */
int find_letter_source(LPOLetter_T *letter,int iseq)
{
  LPOLetterSource_T *source;
  for (source= &letter->source;source;source=source->more)/*SCAN SOURCES*/
    if (source->iseq == iseq) /*MATCH! */
      return source->ipos; /* RETURN ITS SEQUENCE POSITION */
  return INVALID_LETTER_POSITION;
}





/**@memo finds links to iseq, starting from letter[ipos] and proceeding in the 
    specified direction up to max_step links.  Optionally, the search 
    can be constrained to only match node iseq:match_ipos.  
    If iseq is found, the number 
    of links connecting its letter to letter[ipos] is reported.  Otherwise 
    it returns max_step+1 to indicate iseq was not found within the given 
    distance.  If p_ipos or p_letter_pos are non-NULL, it will return
    the ipos of the position found in iseq, or the index of the found 
    letter[], respectively. */
int find_sequence_link(int iseq,
		       int match_ipos,
		       int start_pos,
		       LPOLetter_T letter[],
		       int max_step,
		       int please_go_right,
		       int *p_ipos,
		       int *p_letter_pos)
{
  int nstep,ipos;
  LPOLetterLink_T *link,*link0;
  
  if (please_go_right) /* CHOOSE THE DESIRED DIRECTION */
    link0= &letter[start_pos].right;
  else
    link0= &letter[start_pos].left;

  for (link=link0;link && link->ipos>=0;link=link->more) /*SCAN ALL LINKS*/
    if ((ipos=find_letter_source(letter+link->ipos,iseq))>=0
	&& (match_ipos<0 /* NO CONSTRAINT ON match_ipos */
	    || match_ipos==ipos)) { /* POSITION MATCHES */
      if (p_ipos) /* HAND BACK THE SEQUENCE ipos TO CALLER */
	*p_ipos = ipos;
      if (p_letter_pos) /* HAND BACK THE letter[] INDEX TO CALLER */
	*p_letter_pos = link->ipos;
      return 1; /* FOUND SEQUENCE ONLY ONE LINK AWAY FROM HERE! */
    }

  if (max_step>1) { /* OK TO TRY ANOTHER LAYER OF RECURSION */
    for (link=link0;link && link->ipos>=0;link=link->more) { /*SCAN ALL LINKS*/
      nstep=find_sequence_link(iseq,match_ipos,link->ipos,letter,max_step-1,
			       please_go_right,p_ipos,p_letter_pos);
      if (nstep<max_step) /* FOUND IT! SO RETURN LINK COUNT */
	return nstep+1;
    }
  }

  return max_step+1; /* INDICATE iseq NOT FOUND IN THE SPECIFIED RANGE */
}



/** constructs a mapping from [iseq,ipos] -> iletter and back.  This is 
  saved as an index seq_to_po[ipos]==>iletter and 
  po_to_seq[iletter]==>ipos .
  At the same time, it also constructs the actual sequence of the 
  individual source sequences from the data stored in the LPO. */
void build_seq_to_po_index(LPOSequence_T *seq)
{
  int i,j;
  LPOLetterSource_T *source;

  LOOP (i,seq->nsource_seq) { /* DUMP EXISTING INDEX, ALLOC NEW*/
    FREE(seq->source_seq[i].seq_to_po);
    FREE(seq->source_seq[i].po_to_seq);
    FREE(seq->source_seq[i].sequence);
    CALLOC(seq->source_seq[i].seq_to_po,seq->source_seq[i].length,int);
    CALLOC(seq->source_seq[i].po_to_seq,seq->length,int);
    CALLOC(seq->source_seq[i].sequence,seq->source_seq[i].length+1,char);
    LOOP (j,seq->length) /*DEFAULT: PO LETTERS DON'T MAP TO ANY LETTER IN SEQ*/
      seq->source_seq[i].po_to_seq[j]= INVALID_LETTER_POSITION;
  }

  LOOP (i,seq->length) { /* MAP EVERY LETTER ONTO SOURCE INDEXES */
    for (source= &seq->letter[i].source;source;source=source->more) {
      seq->source_seq[source->iseq].seq_to_po[source->ipos]=i;/*INVERSE*/
      seq->source_seq[source->iseq].po_to_seq[i]=source->ipos;/*MAPPING*/
      seq->source_seq[source->iseq].sequence[source->ipos]=seq->letter[i].letter;
    }
  }
}




/** CREATES A NEW SOURCEINFO ENTRY ON seq->source_seq[], SAVING THE 
 FIELDS PASSED BY THE CALLER */
int save_lpo_source(LPOSequence_T *seq, 
		     char name[],
		     char title[],
		     int length,
		     int weight,
		     int bundle_id,
		    int ndata,
		    LPONumericData_T data[])
{
  int i,j;

  REBUFF(seq->source_seq,seq->nsource_seq,SOURCE_SEQ_BUFFER_CHUNK,
	 LPOSourceInfo_T);
  i=seq->nsource_seq;
  seq->source_seq[i].title=strdup(title? title:"untitled");
  seq->source_seq[i].length=length;
  seq->source_seq[i].weight=weight; /* DEFAULT WEIGHTING */
  seq->source_seq[i].bundle_id= bundle_id;
  STRNCPY(seq->source_seq[i].name,name,SEQUENCE_NAME_MAX);

  LOOPF(j,ndata) /* SAVE NUMERIC DATA FOR THIS SOURCE */
    cp_numeric_data(seq->source_seq+i,data+j);

  return seq->nsource_seq++; /* INCREMENT SOURCE SEQUENCE LIST COUNT */
}



/** SAVE source_seq[] ENTRIES FROM ONE LPO TO ANOTHER */
int *save_lpo_source_list(LPOSequence_T *seq,
			  int nsource_seq,
			  LPOSourceInfo_T source_seq[])
{ 
  int i,*list=NULL;

  CALLOC(list,nsource_seq,int);
  LOOPF (i,nsource_seq)  
    list[i]=save_lpo_source(seq,source_seq[i].name,source_seq[i].title,
			    source_seq[i].length,source_seq[i].weight,
			    source_seq[i].bundle_id,
			    source_seq[i].ndata,source_seq[i].data);

  return list; /*HAND BACK LIST OF NEW INDICES ASSIGNED TO THESE source_seq*/
}



/** ADDS A LINK TO ipos TO THE LINKED LIST STORED IN list,
 ALLOCATING A NEW LPOLetterLink IF NEEDED */
LPOLetterLink_T *add_lpo_link(LPOLetterLink_T *list,LPOLetterRef_T ipos)
{
  if (list->ipos <0) { /* FIRST ENTRY IS EMPTY, SO USE IT */
    list->ipos = ipos;
    return list;/*RETURNS PTR TO LINK IN WHICH ipos STORED */
  }
  do { /* SCAN LIST CHECKING IF ALREADY STORED */
    if (list->ipos == ipos) /* ALREADY STORED, NO NEED TO STORE AGAIN */
      return list;/*RETURNS PTR TO LINK IN WHICH ipos STORED */
  } while (list->more? (list=list->more):0);

  CALLOC(list->more,1,LPOLetterLink_T); /* ADD ENTRY TO LINKED LIST */
  list->more->ipos=ipos; /* SAVE THE LETTER REFERENCE */
  return list->more;/*RETURNS PTR TO LINK IN WHICH ipos STORED */
}




void add_lpo_sources(LPOLetterSource_T *new_s,LPOLetterSource_T *old_s,
		     int iseq_new[])/*TRANSLATION TO NEW source_seq[] INDEX*/
{
  for (;new_s->more;new_s=new_s->more);/* GO TO END*/
  for (;old_s;old_s=old_s->more) {/*SAVE SOURCES*/
    if (new_s->ipos>=0) { /* ALREADY A SOURCE HERE, SO CREATE NEW ENTRY */
      CALLOC(new_s->more,1,LPOLetterSource_T);
      new_s=new_s->more;
    }
    new_s->iseq=iseq_new[old_s->iseq]; /* SAVE SEQUENCE ID, POSITION */
    new_s->ipos=old_s->ipos;
  }
}


void reindex_lpo_source_seqs (LPOSequence_T *seq, int *perm)
{
  int i, j, len, nseq, *map, *invmap;
  LPOSourceInfo_T tmp;
  LPOLetterSource_T *src;
  
  len = seq->length;
  nseq = seq->nsource_seq;
  
  CALLOC (map, len, int);
  CALLOC (invmap, len, int);
  
  /* BUILD INITIAL MAP AND INVERSE MAP: */
  for (i=0; i<nseq; i++) {
    map[i] = invmap[i] = -1;
  }
  for (i=0; i<nseq; i++) {
    map[i] = perm[i];
    IF_GUARD (map[i]>=nseq || map[i]<0 || invmap[map[i]]!=-1, 1.1, (ERRTXT,"Bad argument! 'perm' must be a permutation of [0,%d]\n",nseq-1), CRASH);
    invmap[map[i]] = i;
  }
  
  /* RENUMBER SEQS IN 'source_info' ENTRIES FOR ALL LETTERS: */
  for (i=0; i<len; i++) {
    for (src=&(seq->letter[i].source); src!=NULL && src->ipos>=0; src=src->more) {
      src->iseq = map[src->iseq];
    }
  }

  /* SHUFFLE 'source_seq' ENTRIES, IN-PLACE, TO NEW ORDER: */
  /* (THIS DESTROYS map AND invmap.) */
  for (i=0; i<nseq; i++) if ((j = invmap[i]) != i) {
    /* SWAP POSITIONS i AND j: */
    tmp = seq->source_seq[i];
    seq->source_seq[i] = seq->source_seq[j];
    seq->source_seq[j] = tmp;
    /* UPDATE map AND invmap: */
    map[j] = map[i];
    invmap[map[i]] = j;
    map[i] = invmap[i] = i;
  }
  
  FREE (map);
  FREE (invmap);
}


void copy_lpo_letter(LPOLetter_T *new,LPOLetter_T *old,
		     LPOLetterRef_T old_to_new[],
		     int iseq_new[])
{
  LPOLetterLink_T *link;
  new->letter=old->letter; /* SAVE ITS SEQUENCE LETTER */
  add_lpo_sources(&new->source,&old->source,iseq_new); /* SAVE SOURCES */
  for (link= &old->left;link && link->ipos>=0;link=link->more) /*SAVE left*/
    add_lpo_link(&new->left,old_to_new[link->ipos]);
  for (link= &old->right;link && link->ipos>=0;link=link->more)/*SAVE right*/
    add_lpo_link(&new->right,old_to_new[link->ipos]);
  return;
}



 /** FUSES RING a AND RING b BY CROSSLINK OPERATION */
void crosslink_rings(LPOLetterRef_T a,LPOLetterRef_T b,LPOLetter_T seq[])
{
  LPOLetterRef_T align_ring;

  if (seq[a].ring_id==seq[b].ring_id) /* ALREADY ON SAME RING, DO NOTHING! */
    return;
  else if (seq[a].ring_id<seq[b].ring_id) { /* a<b SO RESET b.ring_id */
    align_ring=b; /* TRAVERSE RING b */
    do seq[align_ring].ring_id=seq[a].ring_id; /* ENFORCE a.ring_id == b.ring_id*/
    while ((align_ring=seq[align_ring].align_ring) != b);
  }
  else { /* a>b SO RESET a.ring_id */
    align_ring=a; /* TRAVERSE RING a */
    do seq[align_ring].ring_id=seq[b].ring_id; /* ENFORCE a.ring_id == b.ring_id*/
    while ((align_ring=seq[align_ring].align_ring) != a);
  }

  align_ring=seq[a].align_ring; /* JUST LIKE A SWAP */
  seq[a].align_ring=seq[b].align_ring; /* CROSSLINK THE TWO RINGS */
  seq[b].align_ring=align_ring;
  return;
}




void copy_old_ring_to_new(LPOLetterRef_T start, /* START OF RING TO COPY*/
			  LPOLetter_T old_lpo[],
			  LPOLetter_T new_lpo[],
			  LPOLetterRef_T old_to_new[]) /* MAPPING */
{
  LPOLetterRef_T ipos,next_pos;
  for (ipos=start;(next_pos=old_lpo[ipos].align_ring) != start;ipos=next_pos)
    crosslink_rings(old_to_new[ipos],old_to_new[next_pos],new_lpo);
}



/* TEMPORARY: THESE CONSTANTS CONTROL SEGMENT FUSION
   TESTING THE FOLLOWING NUMBERS ON DNA ASSEMBLY*/
#define FISSION_BREAK 5 /* LENGTH OF MISMATCH THAT SPLITS INTO SEGMENTS */
#define MINIMUM_FUSION 10 /* MINIMUM #IDENTITIES FOR SEGMENT TO FUSE */
#define FUSION_PERCENT 0.8 /* MINIMUM OVERALL IDENTITY FOR SEGMENT TO FUSE*/

void mark_fusion_segments(int len_x,LPOLetter_T seq_x[],
			  int len_y,LPOLetter_T seq_y[],
			  LPOLetterRef_T y_to_x[],
			  int fission_break,
			  int minimum_fusion_length,
			  float minimum_fusion_identity,
			  char do_fuse[])
{
  int i,i_x,i_y,mismatch_length=0,identity_count=0,fission_break_point= -1;

  LOOP (i_y,len_y) {
    if ((i_x=y_to_x[i_y])>=0 && seq_x[i_x].letter==seq_y[i_y].letter)
      do_fuse[i_y]=1;
  }
#ifdef SOURCE_EXCLUDED
  LOOPF (i_y,len_y) {
    if ((i_x=y_to_x[i_y])<0 /* NOT ALIGNED AT ALL */
	|| seq_x[i_x].letter!=seq_y[i_y].letter /* MISMATCH */
	|| i_y==len_y-1) { /* END OF SEQ, MUST CHECK LAST SEGMENT! */
      if (++mismatch_length>=fission_break /* BREAK POINT */
	  || i_y==len_y-1) { /* END OF SEQ, MUST CHECK LAST SEGMENT! */
	if (identity_count>=minimum_fusion_length /* END OF A FUSION SEGMENT?*/
	    && minimum_fusion_identity
	       *(i_y-mismatch_length-fission_break_point)<= identity_count) {
	  /* YES, MARK PRECEEDING SEGMENT FOR FUSION*/
	  for (i=i_y-fission_break;i>fission_break_point;i--)/*MARK SEGMENT!*/
	    if ((i_x=y_to_x[i])>=0 && seq_x[i_x].letter==seq_y[i].letter)
	      do_fuse[i]=1; /* IDENTITY! MARK THIS POSITION TO BE FUSED */
	}
	fission_break_point=i_y; /* NO FUSION SEGMENT CAN EXTEND PAST HERE */
	identity_count=0; /* RESET IDENTITY COUNTER FOR STARTING NEXT SEGMENT*/
      }
    }
    else { /* PERFECT IDENTITY */
      mismatch_length=0;
      identity_count++;
    }
  }
#endif
}




int reindex_lpo_fusion(int len_x,LPOLetter_T seq_x[],
		       int len_y,LPOLetter_T seq_y[],
		       LPOLetterRef_T x_to_y[],
		       LPOLetterRef_T y_to_x[],
		       LPOLetterRef_T new_x[],
		       LPOLetterRef_T new_y[],
		       int fission_break,
		       int minimum_fusion_length,
		       float minimum_fusion_identity)
{
  int new_len=0;
  LPOLetterRef_T i_x,i_y,i_ring,end_of_ring= -1;
  char *do_fuse=NULL;

  CALLOC(do_fuse,len_y,char); /* MARK POSITIONS TO FUSE seq_y TO seq_x */
  mark_fusion_segments(len_x,seq_x,len_y,seq_y,y_to_x,fission_break,
		       minimum_fusion_length,minimum_fusion_identity,do_fuse);

  for (i_x=i_y=new_len=0;i_x<len_x;i_x++) { /* BUILD new_x, new_y */
    for (i_ring=i_x;i_ring<len_x && seq_x[i_ring].ring_id==seq_x[i_x].ring_id;
	 i_ring++) /* SCAN X RING, AND SEE IF ANYTHING ALIGNS TO Y */
      if (x_to_y[i_ring]>=0) { /* IF SO, INSERT Y NOW TO KEEP X RING TOGETHER*/
	while (i_y<x_to_y[i_ring]) /* CONCATENATE LBEFORE ALIGNED LETTER */
	  new_y[i_y++]=new_len++; /* ADD TO new_lpo */
	break;
      }
	
    if (x_to_y[i_x]>=0   /* ALIGNED, SO FIRST INSERT PRECEEDING FROM seq_y */
	&& i_y<len_y) {
      for (i_ring=seq_y[i_y].align_ring;i_ring!=i_y; /* CIRCLE THE RING*/
	   i_ring=seq_y[i_ring].align_ring) /*FIND END OF i_y ALIGNMENT RING*/
	if (i_ring>end_of_ring) /* FIND MAXIMUM INDEX ON THIS RING */
	  end_of_ring=i_ring; /* RING GUARANTEED TO BE ONE CONTIGUOUS BLOCK*/

      if (do_fuse[i_y]) /* THIS POSITION MEETS OUR FUSION CRITERIA, SO FUSE*/
	new_y[i_y++]=new_len; /* USE SAME INDEX AS WILL BE USED FOR i_x */
      else /* NOT IDENTICAL, SO GIVE IT ITS OWN LETTER */
	new_y[i_y++]=new_len++;  /* ADD TO new_lpo */
    }
    new_x[i_x]=new_len++; /* ADD TO new_lpo */

    while (i_y<=end_of_ring) /* CONCATENATE LETTERS ALIGNED TO i_y */
      new_y[i_y++]=new_len++; /*KEEP ALIGNED LETTERS AS ONE CONTIGUOUS BLOCK!*/
  }
  while (i_y<len_y) { /* CONCATENATE TAIL OF seq_y ONTO new_lpo */
    new_y[i_y++]=new_len++; /* ADD TO new_lpo */
  }

  FREE(do_fuse);
  return new_len;
}


void remap_x_to_new(int nremap_x,/*TRANSLATE INDICES IN remap_x USING new_x*/
		    LPOLetterRef_T remap_x[],/* ARRAY OF INDICES TO TRANSLATE*/
		    int len_x,
		    LPOLetterRef_T new_x[]) /* MAPPING FROM OLD -> NEW */
{
  int i;

  LOOP (i,nremap_x) {
    if (remap_x[i]>=0 && remap_x[i]<len_x) /* KEEP IN BOUNDS */
      remap_x[i]=new_x[remap_x[i]];
    else /* BOOM!! OUT OF RANGE!! */
      abort(); /* DELIBERATELY CRASH TO TRAP THE ERROR! */
  }
}



LPOSequence_T *copy_fuse_lpo_remap(LPOSequence_T *holder_x, /*FUSES TWO LPOs*/
			     LPOSequence_T *holder_y, /*MAKES *NEW* HOLDER */
			     LPOLetterRef_T x_to_y[],
			     LPOLetterRef_T y_to_x[],
			     int nremap_x,
			     LPOLetterRef_T remap_x[])
{
  int i,new_len,*iseq_new=NULL,len_x,len_y;
  LPOLetterRef_T i_x,i_y,*new_x=NULL,*new_y=NULL;
  LPOLetter_T *new_lpo=NULL,*seq_x,*seq_y;
  LPOSequence_T *new_seq=NULL;

  len_x=holder_x->length; /* GET letter ARRAY FROM BOTH x AND y */
  seq_x=holder_x->letter;
  len_y=holder_y->length;
  seq_y=holder_y->letter;

  CALLOC(new_seq,1,LPOSequence_T); /* CREATE A NEW HOLDER */
  CALLOC(new_x,len_x,LPOLetterRef_T);
  CALLOC(new_y,len_y,LPOLetterRef_T);
  new_len=reindex_lpo_fusion(len_x,seq_x,len_y,seq_y,x_to_y,y_to_x,new_x,new_y,
			     FISSION_BREAK,MINIMUM_FUSION,FUSION_PERCENT);
  CALLOC(new_lpo,new_len,LPOLetter_T); /* ALLOCATE NEW LINEARIZED PO */
  LOOP (i,new_len) { /* INITIALIZE ALL LINKS TO INVALID */
    new_lpo[i].left.ipos=new_lpo[i].right.ipos=new_lpo[i].source.ipos
      = INVALID_LETTER_POSITION;
    new_lpo[i].align_ring=new_lpo[i].ring_id=i; /* POINT TO SELF */
  }
  new_seq->length=new_len; /* SAVE NEW LPO ARRAY IN NEW HOLDER */
  new_seq->letter=new_lpo;

  iseq_new=save_lpo_source_list(new_seq,holder_x->nsource_seq,/*COPY x SOURCE*/
				holder_x->source_seq);
  LOOP (i_x,len_x) /* COPY LETTER DATA TO CORRESPONDING LETTERS OF NEW LPO */
    copy_lpo_letter(new_lpo+new_x[i_x],seq_x+i_x,new_x,iseq_new);
  FREE(iseq_new);
  iseq_new=save_lpo_source_list(new_seq,holder_y->nsource_seq,/*COPY y SOURCE*/
				holder_y->source_seq);
  LOOP (i_y,len_y)
    copy_lpo_letter(new_lpo+new_y[i_y],seq_y+i_y,new_y,iseq_new);
  FREE(iseq_new);

  LOOP (i_x,len_x) /* COPY OLD ALIGNMENT RINGS TO THE NEW LPO */
    copy_old_ring_to_new(i_x,seq_x,new_lpo,new_x);
  LOOP (i_y,len_y)
    copy_old_ring_to_new(i_y,seq_y,new_lpo,new_y);

  LOOP (i_x,len_x) /* SAVE ALIGNMENT OF seq_x AND seq_y TO new_lpo.align_ring*/
    if (x_to_y[i_x]>=0) /* seq_x[i_x] ALIGNED TO seq_y, SO SAVE! */
      crosslink_rings(new_x[i_x],new_y[x_to_y[i_x]],new_lpo);

  if (remap_x) /* CONVERT OLD INDEX TABLE TO NEW REFERENCE SYSTEM */
    remap_x_to_new(nremap_x,remap_x,len_x,new_x);
  FREE(new_x); /* DUMP SCRATCH MEMORY AND RETURN */
  FREE(new_y);
  return new_seq; /* HAND BACK THE NEW LPO CONTAINING THE FUSION */
}


/** fuses the two partial orders holder_x and holder_y, based upon the 
 letter_x <--> letter_y mapping specified by x_to_y and y_to_x (which
 must be consistent!)  A new LPO is created to store the result;
 neither holder_x or holder_y are changed */
LPOSequence_T *copy_fuse_lpo(LPOSequence_T *holder_x, /*WRAPPER: NO REMAPPING*/
			LPOSequence_T *holder_y,
			LPOLetterRef_T x_to_y[],
			LPOLetterRef_T y_to_x[])
{
  return copy_fuse_lpo_remap(holder_x,holder_y,x_to_y,y_to_x,0,NULL);
}




LPOSequence_T *copy_lpo(LPOSequence_T *holder_x)
{
  int i;
  LPOSequence_T dummy,*new_copy;
  LPOLetterRef_T *x_to_y=NULL;
  
  memset(&dummy,0,sizeof(dummy)); /* BLANK ALL FIELDS */
  CALLOC(x_to_y,holder_x->length,LPOLetterRef_T);
  LOOP (i,holder_x->length)
    x_to_y[i]= INVALID_LETTER_POSITION;

  new_copy=copy_fuse_lpo(holder_x,&dummy,x_to_y,x_to_y);
  FREE(x_to_y);
  return new_copy;
}





void translate_lpo(int old_len,LPOLetterRef_T old_to_new[],
			 LPOLetter_T seq[])
{
  int i,block_end;
  LPOLetterLink_T *link;
  block_end=old_len;
  LOOPB (i,old_len) { /* TRANSLATE AND SHIFT THE WHOLE LPO */
    seq[i].align_ring=old_to_new[seq[i].align_ring]; /* TRANSLATE INDICES*/
    seq[i].ring_id=old_to_new[seq[i].ring_id];
    for (link= &seq[i].left;link && link->ipos>=0;link=link->more)
      link->ipos = old_to_new[link->ipos];
    for (link= &seq[i].right;link && link->ipos>=0;link=link->more)
      link->ipos = old_to_new[link->ipos];

    if ((0==i || /* HIT SEQ START, SO COPY THE BLOCK NOW!*/
	old_to_new[i] != old_to_new[i-1]+1) /*BLOCK BOUNDARY*/
	&& old_to_new[i] != i) { /* ACTUALLY REQUIRES A SHIFT */
      memmove(seq+old_to_new[i],seq+i,(block_end-i)*sizeof(LPOLetter_T));
      block_end=i; /* RESET TO POINT TO END OF NEXT BLOCK TO COPY */
    } /* BLOCK NOW SHIFTED TO ITS NEW LOCATION */
  }
}



LPOSequence_T *fuse_lpo_remap(LPOSequence_T *holder_x,
			LPOSequence_T *holder_y,
			LPOLetterRef_T x_to_y[],
			LPOLetterRef_T y_to_x[],
			int nremap_x,
			LPOLetterRef_T remap_x[])
{
  int i,new_len,*iseq_new=NULL,len_x,len_y;
  LPOLetterRef_T i_x,i_y,*new_x=NULL,*new_y=NULL;
  LPOLetter_T *new_lpo=NULL,*seq_x,*seq_y;

  len_x=holder_x->length; /* GET letter ARRAY FROM BOTH x AND y */
  seq_x=holder_x->letter;
  len_y=holder_y->length;
  seq_y=holder_y->letter;

  CALLOC(new_x,len_x,LPOLetterRef_T);
  CALLOC(new_y,len_y,LPOLetterRef_T);
  new_len=reindex_lpo_fusion(len_x,seq_x,len_y,seq_y,x_to_y,y_to_x,new_x,new_y,
			     FISSION_BREAK,MINIMUM_FUSION,FUSION_PERCENT);
  REALLOC(seq_x,new_len,LPOLetter_T); /* EXPAND seq_x TO HOLD FUSED LPO */
  translate_lpo(len_x,new_x,seq_x); /* SHIFT ALL THE LETTERS TO NEW LOCATIONS*/
  new_lpo=seq_x; /* THE EXPANDED VERSION OF seq_x IS OUR NEW LPO */
  holder_x->length=new_len; /* SAVE NEW LPO ARRAY IN NEW HOLDER */
  holder_x->letter=new_lpo;

  LOOP (i_y,len_y) { /* INITIALIZE ALL LETTERS FOR STORING seq_y TO BLANK */
    if (y_to_x[i_y]>=0 && new_x[y_to_x[i_y]]==new_y[i_y])/*i_y FUSED TO i_x*/
      continue; /* THIS LETTER IS ALREADY PART OF seq_x SO DON'T OVERWRITE!!*/
    i=new_y[i_y]; /* TRANSLATE TO NEW INDEXING */
    memset(new_lpo+i,0,sizeof(LPOLetter_T)); /* NULL INITIALIZE IT! */
    new_lpo[i].left.ipos=new_lpo[i].right.ipos=new_lpo[i].source.ipos
      = INVALID_LETTER_POSITION; /* RESET TO UNLINKED STATE */
    new_lpo[i].align_ring=new_lpo[i].ring_id=i; /* POINT TO SELF */
  }

  iseq_new=save_lpo_source_list(holder_x,holder_y->nsource_seq,/*COPY y SRC*/
				holder_y->source_seq);
  LOOP (i_y,len_y) /* COPY LETTER DATA TO CORRESPONDING LETTERS OF NEW LPO */
    copy_lpo_letter(new_lpo+new_y[i_y],seq_y+i_y,new_y,iseq_new);
  FREE(iseq_new);

  LOOP (i_y,len_y) /* COPY OLD ALIGNMENT RINGS TO THE NEW LPO */
    copy_old_ring_to_new(i_y,seq_y,new_lpo,new_y);

  LOOP (i_x,len_x) /* SAVE ALIGNMENT OF seq_x AND seq_y TO new_lpo.align_ring*/
    if (x_to_y[i_x]>=0) /* seq_x[i_x] ALIGNED TO seq_y, SO SAVE! */
      crosslink_rings(new_x[i_x],new_y[x_to_y[i_x]],new_lpo);

  if (remap_x) /* CONVERT OLD INDEX TABLE TO NEW REFERENCE SYSTEM */
    remap_x_to_new(nremap_x,remap_x,len_x,new_x);
  FREE(new_x); /* DUMP SCRATCH MEMORY AND RETURN */
  FREE(new_y);
  return holder_x; /* HAND BACK x LPO CONTAINING THE FUSION */
}


/** fuses the two partial orders holder_x and holder_y, based upon the 
 letter_x <--> letter_y mapping specified by x_to_y and y_to_x (which
 must be consistent!)  The result is returned in holder_x */
LPOSequence_T *fuse_lpo(LPOSequence_T *holder_x, /*WRAPPER: NO REMAPPING*/
			LPOSequence_T *holder_y,
			LPOLetterRef_T x_to_y[],
			LPOLetterRef_T y_to_x[])
{
  return fuse_lpo_remap(holder_x,holder_y,x_to_y,y_to_x,0,NULL);
}



/** FREES the linked list link including all nodes beneath it; NB: link 
 itself is freed, so DO NOT pass a static LPOLetterLink */
void free_lpo_link_list(LPOLetterLink_T *link)
{
  LPOLetterLink_T *next;
  for (;link;link=next) { /* DUMP ALL THE LINKS */
    next=link->more;
    free(link);
  }
}

/** FREES the linked list source including all nodes beneath it; NB: source 
 itself is freed, so DO NOT pass a static LPOLetterSource */
void free_lpo_source_list(LPOLetterSource_T *source)
{
  LPOLetterSource_T *next;
  for (;source;source=next) { /* DUMP ALL THE SOURCE ENTRIES */
    next=source->more;
    free(source);
  }
}


/** FREES ALL DATA ASSOCIATED WITH letter[], AND OPTIONALLY letter ITSELF*/
void free_lpo_letters(int nletter,LPOLetter_T *letter,int please_free_block)
{
  int i;
  if (!letter) /* NOTHING TO FREE... */
    return;
  LOOP (i,nletter) { /* DUMP ALL LINKED LISTS */
    if (letter[i].left.more)
      free_lpo_link_list(letter[i].left.more);
    if (letter[i].right.more)
      free_lpo_link_list(letter[i].right.more);
    if (letter[i].source.more)
      free_lpo_source_list(letter[i].source.more);
  }
  if (please_free_block) /*DON'T ALWAYS WANT TO FREE... MIGHT BE IN AN ARRAY*/
    free(letter);
}




void free_lpo_sourceinfo(int nsource_seq,LPOSourceInfo_T *source_seq,
			int please_free_block)
{
  int i;
  LOOP (i,nsource_seq) {
    FREE(source_seq[i].title);
    FREE(source_seq[i].sequence);
    FREE(source_seq[i].seq_to_po);
    FREE(source_seq[i].po_to_seq);
    free_lpo_numeric_data(source_seq[i].ndata,source_seq[i].data,TRUE);
    source_seq[i].data=NULL; /* DON'T LEAVE DANGLING POINTER! */
  }
  if (please_free_block)
    free(source_seq);
}




/** FREES ALL DATA FROM seq, AND OPTIONALLY seq ITSELF */
void free_lpo_sequence(LPOSequence_T *seq,int please_free_holder)
{
  int i;
  if (!seq) /* NOTHING TO FREE... */
    return;
  free_lpo_letters(seq->length,seq->letter,TRUE);
  seq->letter=NULL; /* MARK AS FREED... DON'T LEAVE DANGLING POINTER! */
  FREE(seq->title);
  FREE(seq->sequence);
  if (seq->source_seq) {
    free_lpo_sourceinfo(seq->nsource_seq,seq->source_seq,TRUE);
    seq->source_seq=NULL; /* MARK AS FREED... DON'T LEAVE DANGLING POINTER! */
  }
  if (please_free_holder) /*DON'T ALWAYS WANT TO FREE... MIGHT BE IN AN ARRAY*/
    free(seq);
}

/**@memo {\bfEXAMPLE}: dump the LPO dna_lpo, and all its associated data: \begin{verbatim}
  if (dna_lpo) 
    free_lpo_sequence(dna_lpo,TRUE);
\end{verbatim}
*/


/** creates a new sequence which follows path[] through the LPO seq,
 and gives it the specified name and title  */
int add_path_sequence(int path_length,
		      LPOLetterRef_T path[],
		      LPOSequence_T *seq,
		      char name[],
		      char title[])
{
  int i,iseq_new;
  LPOSourceInfo_T *new_seq;
  LPOLetterSource_T save_source={0,0,NULL};

  /* CREATE SOURCE ENTRY FOR THIS SEQUENCE */
  iseq_new=save_lpo_source(seq,name,title,path_length,0,NO_BUNDLE,0,NULL);

  LOOP (i,path_length) { /* ADD THIS AS SOURCE TO ALL POSITIONS IN path */
    save_source.ipos=i;
    add_lpo_sources(&seq->letter[path[i]].source,&save_source,&iseq_new);
  } /* NB: THIS DOESN'T CHECK THAT path IS A VALID WALK THRU THE PARTIAL ORDER
       MIGHT BE A GOOD IDEA TO CATCH POSSIBLE ERRORS IN path */
  return iseq_new; /* RETURN INDEX OF NEWLY CREATED ENTRY */
}

/**@memo EXAMPLE: create a consensus sequence from a path: 
    iseq=add_path_sequence(path_length,path,seq,name,title); 
------------------------------------------------------- 
-------------------------------------------------
*/


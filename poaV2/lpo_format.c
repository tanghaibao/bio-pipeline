

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/** writes the LPO in seq to the stream ifile; optionally a symbol table
 may be given for translating the letters in the LPO to text */
void write_lpo(FILE *ifile,LPOSequence_T *seq,
	       ResidueScoreMatrix_T *score_matrix)
{
  int i;
  LPOLetterLink_T *link;
  LPOLetterSource_T *source;

  fprintf(ifile,"VERSION=LPO.0.1\n");
  fprintf(ifile,"NAME=%s\nTITLE=%s\nLENGTH=%d\nSOURCECOUNT=%d\n",
	  seq->name,seq->title,seq->length,seq->nsource_seq);

  LOOPF (i,seq->nsource_seq)
    fprintf(ifile,"SOURCENAME=%s\nSOURCEINFO=%d %d %d %d %s\n",
	    seq->source_seq[i].name,seq->source_seq[i].length,
	    seq->source_seq[i].istart,seq->source_seq[i].weight,
	    seq->source_seq[i].bundle_id,seq->source_seq[i].title);

  LOOPF (i,seq->length) {
    fprintf(ifile,"%c:",
	    seq->letter[i].letter < score_matrix->nsymbol ? 
	    score_matrix->symbol[seq->letter[i].letter] 
	    : seq->letter[i].letter);
    for (link= &seq->letter[i].left;link && link->ipos>=0;link=link->more)
      fprintf(ifile,"L%d",link->ipos);
    for (source= &seq->letter[i].source;source;source=source->more)
      fprintf(ifile,"S%d",source->iseq); /* SOURCE ID */
    if (seq->letter[i].align_ring!=i) /* ALIGNED TO SOMETHING ELSE */
      fprintf(ifile,"A%d",seq->letter[i].align_ring);
    fputc('\n',ifile);
  }
}
/**@memo example: writing a PO file: 
  if (lpo_file_out) 
    write_lpo(lpo_file_out,lpo_out,score_matrix.symbol);
*/




/** reads an LPO from the stream ifile, dynamically allocates memory for
it, and returns a pointer to the LPO */ 
LPOSequence_T *read_lpo(FILE *ifile)
{
  int i,j,length,nsource_seq,istart,field_id,*pos_count=NULL,value;
  int weight,bundle_id,last_alloc=0;
  LPOSequence_T *seq=NULL;
  char c,name[1024]="",title[4096]="",version[256]="";
  LPOLetterSource_T save_source={0,0,NULL};

  CALLOC(seq,1,LPOSequence_T);
  fscanf(ifile,"VERSION=%s",version);
  fscanf(ifile," NAME=%[^\n]",name);
  fscanf(ifile," TITLE=%[^\n]",title);
  if (fscanf(ifile," LENGTH=%d SOURCECOUNT=%d",
	     &length,&nsource_seq)!=2)
    return NULL;
  STRNCPY(seq->name,name,SEQUENCE_NAME_MAX);
  seq->title=strdup(title);
  seq->length=length;
  seq->nsource_seq=nsource_seq;
  CALLOC(seq->letter,length,LPOLetter_T);
  GETMEM(seq->source_seq,nsource_seq,last_alloc,SOURCE_SEQ_BUFFER_CHUNK,LPOSourceInfo_T);
  CALLOC(pos_count,nsource_seq,int);
  LOOP (i,length) { /* INITIALIZE ALL LINKS TO INVALID */
    seq->letter[i].align_ring=i; /* POINT TO SELF */
    seq->letter[i].ring_id= INVALID_LETTER_POSITION; /* BLANK! */
    seq->letter[i].left.ipos=seq->letter[i].right.ipos=
      seq->letter[i].source.ipos= INVALID_LETTER_POSITION;
  }

  LOOPF(i,nsource_seq) { /* SAVE SOURCE INFO LIST */
    if (fscanf(ifile," SOURCENAME=%[^\n] SOURCEINFO=%d %d %d %d",
	       name,&length,&istart,&weight,&bundle_id)!=5)
      return NULL;
    /* SKIP WHITESPACE BEFORE TITLE; ALLOW EMPTY TITLE. */
    fscanf(ifile,"%*[ \t]");
    if (fscanf(ifile,"%[^\n]",title)!=1)
      title[0] = '\0';
    STRNCPY(seq->source_seq[i].name,name,SEQUENCE_NAME_MAX);
    seq->source_seq[i].length=length;
    seq->source_seq[i].istart=istart;
    seq->source_seq[i].weight=weight;
    seq->source_seq[i].bundle_id=bundle_id;
    seq->source_seq[i].title=strdup(title);
  }

  LOOPF (i,seq->length) { /* NOW READ THE ACTUAL PARTIAL ORDER */
    if (fscanf(ifile," %c:",&c)!=1) /* READ SEQUENCE LETTER */
      return NULL;
    seq->letter[i].letter=c;
    while ((field_id=getc(ifile))!=EOF && '\n'!=field_id) {/* READ FIELDS*/
      if (1!=fscanf(ifile,"%d",&value))
	return NULL;
      switch (field_id) {
      case 'L':
	add_lpo_link(&seq->letter[i].left,value); /* ADD LEFT-RIGHT LINKS*/
	add_lpo_link(&seq->letter[value].right,i);
	break;
      case 'S':  /* SAVE THE SOURCE ID */
	save_source.ipos=pos_count[value]++;
	add_lpo_sources(&seq->letter[i].source,&save_source,&value);
	break;
      case 'A': /* SAVE THE ALIGN RING POINTER */
	seq->letter[i].align_ring=value;
	break;
      }
    }
  }

  LOOPF (i,seq->length) { /* SET ring_id TO MINIMUM VALUE ON EACH RING */
    if (seq->letter[i].ring_id<0) {/* NEW RING, UPDATE IT! */
      j=i; /* GO AROUND THE ENTIRE RING, SETTING ring_id TO i */
      do seq->letter[j].ring_id=i; /* i IS MINIMUM VALUE ON THIS RING */
      while ((j=seq->letter[j].align_ring)!=i);
    }
  }

  FREE(pos_count);
  return seq;
}



#define INVALID_LPO_LINK (-99)

enum {
  default_retention_mode,
  default_no_retention_mode
};

/** reads an LPO from the stream ifile, dynamically allocates memory for
it, and returns a pointer to the LPO */ 
LPOSequence_T *read_lpo_select(FILE *ifile,FILE *select_file,
			       int keep_all_links,int remove_listed_sequences)
{
  int i,j,k,length,nsource_seq,istart,field_id,*pos_count=NULL,value;
  int weight,bundle_id,last_alloc=0,*iseq_compact=NULL,*last_pos=NULL;
  int nlink,*link_list=NULL,*match_pos=NULL,*ring_old=NULL;
  int *pos_compact=NULL,npos_compact=0,keep_this_letter,retention_mode;
  LPOSequence_T *seq=NULL;
  char c,name[1024]="",title[4096]="",version[256]="";
  LPOLetterSource_T save_source={0,0,NULL},*source=NULL;

  if (remove_listed_sequences)
    retention_mode=default_retention_mode;/*KEEP SEQS AS DFLT, SKIP IF LISTED*/
  else /* SKIP SEQS UNLESS LISTED IN select_file */
    retention_mode=default_no_retention_mode;

  CALLOC(seq,1,LPOSequence_T);
  fscanf(ifile,"VERSION=%s",version);
  fscanf(ifile," NAME=%[^\n]",name);
  fscanf(ifile," TITLE=%[^\n]",title);
  if (fscanf(ifile," LENGTH=%d SOURCECOUNT=%d",
	     &length,&nsource_seq)!=2)
    return NULL;
  STRNCPY(seq->name,name,SEQUENCE_NAME_MAX);
  seq->title=strdup(title);
  seq->length=length;
  seq->nsource_seq=nsource_seq;
  CALLOC(seq->letter,length,LPOLetter_T);
  GETMEM(seq->source_seq,nsource_seq,last_alloc,SOURCE_SEQ_BUFFER_CHUNK,LPOSourceInfo_T);
  CALLOC(pos_count,nsource_seq,int);
  LOOP (i,length) { /* INITIALIZE ALL LINKS TO INVALID */
    seq->letter[i].align_ring=i; /* POINT TO SELF */
    seq->letter[i].ring_id= INVALID_LETTER_POSITION; /* BLANK! */
    seq->letter[i].left.ipos=seq->letter[i].right.ipos=
      seq->letter[i].source.ipos= INVALID_LETTER_POSITION;
  }

  LOOPF(i,nsource_seq) { /* SAVE SOURCE INFO LIST */
    if (fscanf(ifile," SOURCENAME=%[^\n] SOURCEINFO=%d %d %d %d",
	       name,&length,&istart,&weight,&bundle_id)!=5)
      return NULL;
    /* SKIP WHITESPACE BEFORE TITLE; ALLOW EMPTY TITLE. */
    fscanf(ifile,"%*[ \t]");
    if (fscanf(ifile,"%[^\n]",title)!=1)
      title[0] = '\0';
    STRNCPY(seq->source_seq[i].name,name,SEQUENCE_NAME_MAX);
    seq->source_seq[i].length=length;
    seq->source_seq[i].istart=istart;
    seq->source_seq[i].weight=weight;
    seq->source_seq[i].bundle_id=bundle_id;
    seq->source_seq[i].title=strdup(title);
  }

  CALLOC(iseq_compact,nsource_seq,int);
  CALLOC(last_pos,nsource_seq,int);
  CALLOC(match_pos,nsource_seq,int);
  if (select_file) {
    LOOP (i,nsource_seq) /* DEFAULT: MARKED AS INVALID */
      iseq_compact[i]= -retention_mode;
    while (fscanf(select_file,"SOURCENAME=%[^\n]\n",name)==1) {
      LOOP (i,nsource_seq)
	if (strcmp(seq->source_seq[i].name,name)==0) {
	  iseq_compact[i]= retention_mode-1;
	  break;
	}
    }
  }

  j=0;
  LOOPF (i,nsource_seq) {
    last_pos[i]= -1; /* DEFAULT: INVALID */
    if (iseq_compact[i]>=0) {
      iseq_compact[i]=j;
      if (i>j)
	memcpy(seq->source_seq+j,seq->source_seq+i,sizeof(LPOSourceInfo_T));
      match_pos[j]= INVALID_LPO_LINK; /* DEFAULT: NO VALID LINK! */
      j++;
    }
    else { /* EXCLUDE THIS SEQUENCE FROM THE FILTERED LPO */
      iseq_compact[i]= INVALID_LETTER_POSITION;
      if (seq->source_seq[i].title) 
	free(seq->source_seq[i].title);
    }
  }
  seq->nsource_seq=nsource_seq=j; /* COMPACTED COUNT OF SEQUENCES TO KEEP*/

  CALLOC(link_list,seq->length,int); /*TEMPORARY DATA FOR COMPACTION MAPPING */
  CALLOC(pos_compact,seq->length,int);
  CALLOC(ring_old,seq->length,int);
  npos_compact=0;

  LOOPF (i,seq->length) { /* NOW READ THE ACTUAL PARTIAL ORDER */
    if (fscanf(ifile," %c:",&c)!=1) /* READ SEQUENCE LETTER */
      return NULL;
    seq->letter[npos_compact].letter=c;
    nlink=0;
    keep_this_letter=0; /*DEFAULT */
    while ((field_id=getc(ifile))!=EOF && '\n'!=field_id) {/* READ FIELDS*/
      if (1!=fscanf(ifile,"%d",&value))
	return NULL;
      switch (field_id) {
      case 'L':
	if (pos_compact[value]>=0) /*COULD BE VALID LINK: WAIT TO CHECK SRCs*/
	  link_list[nlink++]=pos_compact[value]; /* SAVE IT TEMPORARILY */
	break;
      case 'S':  /* SAVE THE SOURCE ID */
	if (iseq_compact[value]>=0) { /*KEEP THIS SOURCE, SO KEEP THIS LETTER*/
	  keep_this_letter=1;
	  value=iseq_compact[value]; /* TRANSLATE TO ITS COMPACTED INDEX*/
	  if (last_pos[value]>=0) { /* MAKE SURE WE HAVE LINK TO LAST POSITION */
	    LOOP (j,nlink) /* CHECK TO SEE IF LINK ALREADY SAVED */
	      if (link_list[j]==last_pos[value])
		break;
	    if (LOOP_FINISHED(j,nlink)) { /* NO LINK??? ADD IT!!! */
	      link_list[nlink++]=last_pos[value];
	    }
	  }
	  last_pos[value]=npos_compact; /* THIS SEQ POS IS AT THIS NODE */
	  save_source.ipos=pos_count[value]++;/* COUNT LENGTH OF THIS SEQ */
	  add_lpo_sources(&seq->letter[npos_compact].source,
			  &save_source,&value);
	  seq->letter[npos_compact].ring_id= INVALID_LETTER_POSITION;
	  seq->letter[npos_compact].align_ring=i;
	  ring_old[i]=i; /* DEFAULT: SELF-RING OF ONE LETTER*/
	}
	break;
      case 'A': /* SAVE THE ALIGN RING POINTER */
	ring_old[i]=value; /* SAVE OLD ALIGN RING INDICES */
	if (keep_this_letter) /* TEMP'Y: SAVE REVERSE MAPPING TO A-R INDICES*/
	  seq->letter[npos_compact].align_ring=i;
	break;
      }
    }
    if (keep_this_letter) {
      for (source= &seq->letter[npos_compact].source;source;source=source->more)
	match_pos[source->iseq]=source->ipos - 1; /*VALID LINK MUST MATCH m_p*/

      LOOPF (j,nlink) { /*ADD LEFT-RIGHT LINKS*/
	for (source= &seq->letter[link_list[j]].source;source;source=source->more)
	  if (keep_all_links /* KEEP LINKS EVEN IF NOT FROM SELECTED SEQS*/
	      || source->ipos == match_pos[source->iseq]) { /*VALID LINK! */
	    add_lpo_link(&seq->letter[npos_compact].left,link_list[j]);
	    add_lpo_link(&seq->letter[link_list[j]].right,npos_compact);
	    break; /* SAVED THIS LINK! */
	  }
      }
      for (source= &seq->letter[npos_compact].source;source;source=source->more)
	match_pos[source->iseq]= INVALID_LPO_LINK;/*DFLT:NO VALID LINK*/

      pos_compact[i]=npos_compact++;
    }
    else 
      pos_compact[i]= INVALID_LETTER_POSITION;
  }
  seq->length=npos_compact;

  LOOPF (i,seq->length) { /* SET ring_id TO MINIMUM VALUE ON EACH RING */
    if (seq->letter[i].ring_id<0) {/* NEW RING, UPDATE IT! */
      j=seq->letter[i].align_ring; /* GO AROUND THE ENTIRE RING, SETTING ring_id TO i */
      do {
	/*printf("i=%d\tj=%d\tpos_compact[j]=%d\tring_old[j]=%d\n",i,j,pos_compact[j],ring_old[j]);*/
	if (pos_compact[j]>=0) {
	  k=pos_compact[j];
	  seq->letter[k].ring_id=i; /* i IS MINIMUM VALUE ON THIS RING */
	}
	j=ring_old[j]; /* ADVANCE TO NEXT LETTER ON THE RING */
	if (pos_compact[j]>=0) /* IF VALID, POINT IT BACK TO PREVIOUS LETTER*/
	  seq->letter[pos_compact[j]].align_ring=k;
      }
      while (pos_compact[j]!=i); /* STOP WHEN WE'VE COMPLETED THE RING, BACK TO START*/
    }
  }

  FREE(pos_count);
  FREE(pos_compact);
  FREE(iseq_compact);
  FREE(last_pos);
  FREE(link_list);
  FREE(ring_old);
  FREE(match_pos);
  return seq;
}

/*
LPOSequence_T *read_lpo(FILE *ifile)
{
  return read_lpo_select(ifile,NULL);
}
*/




#define FASTA_GAP_CHARACTER '.'
int xlate_lpo_to_al(LPOSequence_T *seq,
		    int nsymbol,char symbol[],int ibundle,
		    char gap_character,
		    char ***p_seq_pos,char **p_p,char **p_include)
{
  int i,j,iring=0,nring=0,current_ring=0,iprint;
  char **seq_pos=NULL,*p=NULL,*include_in_save=NULL;
  LPOLetterSource_T *source;

  LOOPF (i,seq->length) /* COUNT TOTAL #ALIGNMENT RINGS IN THE LPO */
    if (seq->letter[i].ring_id != current_ring) { /* NEXT RING */
      current_ring=seq->letter[i].ring_id;
      nring++;
    }
  nring++; /* DON'T FORGET TO COUNT THE LAST RING!!! */
  
  CALLOC(seq_pos,seq->nsource_seq,char *); /* ALLOCATE MAP ARRAY*/
  CALLOC(p,seq->nsource_seq*nring,char);
  LOOP (i,seq->nsource_seq) /* BUILD POINTER ARRAY INTO MAP ARRAY */
    seq_pos[i]=p+i*nring;
  memset(p,gap_character,seq->nsource_seq*nring);
  /* DEFAULT IS NO SEQUENCE PRESENT AT THIS POSITION */

  current_ring=0; /* RESET TO BEGINNING */
  LOOPF (i,seq->length) { /* NOW MAP THE LPO TO A FLAT LINEAR ORDER */
    if (seq->letter[i].ring_id != current_ring) { /* NEXT RING */
      current_ring=seq->letter[i].ring_id;
      iring++;
    } /* MAP EACH SOURCE SEQ ONTO LINEAR ORDER INDEXED BY iring */
    for (source= &seq->letter[i].source;source;source=source->more)
      if (symbol && seq->letter[i].letter<nsymbol)  /* TRANSLATE TO symbol */
	seq_pos[source->iseq][iring]= symbol[seq->letter[i].letter];
      else  /* NO NEED TO TRANSLATE */
	seq_pos[source->iseq][iring]= seq->letter[i].letter;
  }

  if (ibundle>=0) { /* ONLY SAVE SEQS THAT ARE IN THIS BUNDLE */
    CALLOC(include_in_save,nring,char); /* BLANK FLAGS: WHAT RINGS TO SHOW*/
    LOOP (iring,nring) { /* CHECK EACH RING TO SEE IF IT'S IN BUNDLE */
      LOOP (i,seq->nsource_seq) {
	if (seq_pos[i][iring]!=gap_character /* ALIGNED HERE! */
	    && seq->source_seq[i].bundle_id == ibundle) { /* PART OF BUNDLE!*/
	  include_in_save[iring]=1; /* SO INCLUDE THIS RING */
	  break;
	}
      }
    }
  }
  if (p_seq_pos)
    *p_seq_pos = seq_pos;
  return nring;
}


/** writes the LPO in FASTA format, including all sequences in the 
  specified bundle */
void write_lpo_bundle_as_fasta(FILE *ifile,LPOSequence_T *seq,
			       int nsymbol,char symbol[],int ibundle)
{
  int i,j,nring=0,iprint;
  char **seq_pos=NULL,*p=NULL,*include_in_save=NULL;

  nring=xlate_lpo_to_al(seq,nsymbol,symbol,ibundle, /* TRANSLATE TO */
			FASTA_GAP_CHARACTER,        /* RC-MSA FMT */
			&seq_pos,&p,&include_in_save);
  LOOPF (i,seq->nsource_seq) { /* NOW WRITE OUT FASTA FORMAT */
    if (ibundle<0 /* PRINT ALL BUNDLES */
	|| seq->source_seq[i].bundle_id == ibundle) { /* OR JUST THIS BUNDLE*/
      fprintf(ifile,">%s %s",seq->source_seq[i].name,seq->source_seq[i].title);
      iprint=0;
      LOOPF (j,nring) { /* WRITE OUT 60 CHARACTER SEQUENCE LINES */
	if (NULL==include_in_save || include_in_save[j]) {
	  fprintf(ifile,"%s%c",iprint%60? "":"\n", seq_pos[i][j]);
	  iprint++; /* KEEP COUNT OF PRINTED CHARACTERS */
	}
      }
      fputc('\n',ifile);
    }
  }

  FREE(p); /* DUMP TEMPORARY MEMORY */
  FREE(include_in_save);
  FREE(seq_pos);
}
/**@memo example: writing FASTA format file: 
    if (seq_ifile=fopen(fasta_out,"w")) {
      write_lpo_bundle_as_fasta(seq_ifile,lpo_out,
      score_matrix.nsymbol,score_matrix.symbol,ibundle);
      fclose(seq_ifile);
    }
*/


/** writes the LPO in FASTA format, including all sequences in all bundles 
-------------------------------------------------------
---------------------------------------------------------------------------
*/
void write_lpo_as_fasta(FILE *ifile,LPOSequence_T *seq,
			int nsymbol,char symbol[])
{ /* WRAPPER FUNCTION FOR SAVING ALL BUNDLES!! */
  write_lpo_bundle_as_fasta(ifile,seq,nsymbol,symbol,ALL_BUNDLES);
}




/****************************************************************
  *
  *     WRITE_SEQUENCES
  *
  *     This function writes out sequences.  Each line of sequnces
  *     will be written out in blocks with spacing in between
  *     each block.  
  *
  ***************************************************************/

int write_sequences(FILE *ifile,
		    LPOSequence_T *seq,
		    int indent,
		    int nblock,  /* # OF BLOCKS */
		    int block_size, /* SIZE OF BLOCKS */
		    int block_spacing, /* SPACING BETWEEN BLOCKS */
		    int paragraph_spacing, /* SPACING BETWEEN PRAGRAPHS */
		    int names, /* BOOLEAN: PRINT NAMES AFTER 1ST PARAGRAPH? */
		    char gap_char,
		    int nsymbol,char symbol[])
{
  int i,ip,ial,iseq,iblock,nparagraph,remainder,ipos,ipos_local,broken;
  int len,nring;
  char **seq_pos=NULL,*p=NULL,*include_in_save=NULL;

  nring=xlate_lpo_to_al(seq,nsymbol,symbol,ALL_BUNDLES, /* TRANSLATE TO */
		        gap_char,                   /* RC-MSA FMT */
			&seq_pos,&p,&include_in_save);

  nparagraph = nring/(nblock*block_size);/*# OF FULL PARAGRAPHS */
  remainder = nring%(nblock*block_size);
  if (remainder != 0)
    nparagraph++;
  LOOPF(ip,nparagraph){  /* FOR EACH PARAGRAPH */
    LOOPF(iseq,seq->nsource_seq){ /* FOR EACH SEQUENCE */
      if (ip == 0 || names){
	if ((len=strlen(seq->source_seq[iseq].name))>indent-1){ /*MUST WE TRUNCATE NAME? */
	  LOOPF(i,indent-1) /* PRINT AS MUCH OF NAME AS POSSIBLE  */
	    putc(seq->source_seq[iseq].name[i],ifile);
	  putc(' ',ifile); /* LEAVE A SPACE */
	}
	else{ /* PRINT WHOLE NAME */
	  fprintf(ifile,"%s",seq->source_seq[iseq].name); /* PRINT NAME */
	  LOOP(i,indent-len) /* INDENT LINE */
	    putc(' ',ifile);
	}
      }
      else{
	LOOP(i,indent) /* INDENT LINE */
	  putc(' ',ifile);
      }
      LOOPF(iblock,nblock){ /* FOR EACH BLOCK */
	broken=0;
	LOOPF(i,block_size){
	  ipos = i + (iblock*block_size) + (ip*nblock*block_size);
	  if (ipos>=nring){
	    broken=1;
	    break;
	  }
	  putc(seq_pos[iseq][ipos],ifile); /* APPROP SYMBOL FOR THIS POS*/
	}
	if (broken)
	  break;
	LOOP(i,block_spacing) /* ADD SPACING BETWEEN BLOCKS */
	  putc(' ',ifile);
      } /* END OF BLOCK LOOP */
      putc('\n',ifile); /* NEW LINE AT END OF SEQUENCE LINE */
    } /* END OF SEQ LOOP */
    LOOP(i,paragraph_spacing) /* ADD SPACING BETWEEN PARAGRAPHS */
      putc('\n',ifile);
  } /* END OF PARAGRAPH LOOP */
 done:
  FREE(p); /* DUMP TEMPORARY MEMORY */
  FREE(include_in_save);
  FREE(seq_pos);
  return 0;
}



void export_clustal_seqal(FILE *ifile,
			  LPOSequence_T *seq,
			  int nsymbol,char symbol[])
{
  fprintf(ifile,"CLUSTAL W (1.74) multiple sequence alignment\n\n\n");
  /* WRITE OUT SEQUNENCES: INDENT 36, 1 BLOCK OF 50 CHARS 0 CHARS BLOCK
     SPACING, 2 LINES PARAGRAPH SPACING, PRINT NAMES ON ALL LINES. */
  write_sequences(ifile,seq,36,1,50,0,2,1,'-',nsymbol,symbol);
}


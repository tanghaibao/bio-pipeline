

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"




/** finds the heaviest traversal of the LPO seq[], using dynamic programming;
  at each node the heaviest link is chosen to buildup traversals; finally,
  the traversal with the heaviest overall link weight is returned as an
  array of position indices.  The length of the array is stored in 
  *p_best_len*/
LPOLetterRef_T *heaviest_bundle(int len,LPOLetter_T seq[],
				int nsource_seq,LPOSourceInfo_T source_seq[],
				int *p_best_len)
{
  int i,j,best_right,iright,ibest= -1,best_len=0;
  LPOLetterRef_T *best_path=NULL,*path=NULL;
  LPOLetterLink_T *right;
  LPOLetterSource_T *source;
  LPOScore_T *score=NULL,best_score= -999999,right_score;
  int *contains_pos=NULL,my_overlap,right_overlap;

  CALLOC(path,len,LPOLetterRef_T); /* GET MEMORY FOR DYNAMIC PROGRAMMING */
  CALLOC(score,len,LPOScore_T);
  CALLOC(contains_pos,nsource_seq,int);
  LOOP (i,nsource_seq) /* RESET TO NOT MATCH ANY POSITIONS */
    contains_pos[i]= INVALID_LETTER_POSITION;

  LOOPB (i,len) { /* FIND HEAVIEST PATH BY DYNAMIC PROGRAMMING */
    source= &seq[i].source;  /* MARK SEQUENCES CONTAINING THIS POSITION */
    memset(contains_pos,0,nsource_seq*sizeof(int)); /* ERASE ARRAY */
    do
      if (source_seq[source->iseq].weight>0) /* EXCLUDE ZERO WEIGHT SEQS */
	contains_pos[source->iseq]=source->ipos+1; /* right MUST BE ADJACENT*/
    while (source=source->more); /* KEEP COUNTING TILL NO more */

    right_score=right_overlap=0;  /*DEFAULT MOVE: NOTHING TO THE RIGHT*/
    best_right= INVALID_LETTER_POSITION;
    for (right= &seq[i].right;right && right->ipos>=0;right=right->more) {
      my_overlap=0; /* OVERLAP CALCULATION */
      source= &seq[right->ipos].source;/*COUNT SEQS SHARED IN i AND right*/
      do /* BIAS OVERLAP CALCULATION BY SEQUENCE WEIGHTING */
	if (contains_pos[source->iseq]==source->ipos) /* YES, ADJACENT! */
	  my_overlap += source_seq[source->iseq].weight;
      while (source=source->more); /* KEEP COUNTING TILL NO more */
      
      if (my_overlap>right_overlap /* FIND BEST RIGHT MOVE: BEST OVERLAP */
	  || (my_overlap==right_overlap && score[right->ipos]>right_score)) {
	right_overlap=my_overlap;
	right_score=score[right->ipos];
	best_right=right->ipos;
      }
    }

    path[i]=best_right; /* SAVE THE BEST PATH FOUND */
    score[i]=right_score+right_overlap; /* SAVE THE SCORE */
    if (score[i]>best_score) { /* RECORD BEST SCORE IN WHOLE LPO */
      ibest=i;
      best_score=score[i];
    }
  }

  CALLOC(best_path,len,LPOLetterRef_T); /* MEMORY FOR STORING BEST PATH */
  for (;ibest>=0;ibest=path[ibest])  /* BACK TRACK THE BEST PATH */
    best_path[best_len++]=ibest;

  FREE(path); /* DUMP SCRATCH MEMORY */
  FREE(score);
  FREE(contains_pos);

  if (p_best_len) /* RETURN best_path AND ITS LENGTH */
    *p_best_len = best_len;
  return best_path;
}




int assign_sequence_bundle_id(int path_length,LPOLetterRef_T path[],
			      LPOSequence_T *seq,int bundle_id,
			      float minimum_fraction)
{
  int i,*bundle_count=NULL,nseq_in_bundle=0;
  LPOLetterSource_T *source;
  
  CALLOC(bundle_count,seq->nsource_seq,int);
  LOOP (i,path_length) /* COUNT #POSITIONS OF EACH SEQ ARE IN path */
    for (source= &seq->letter[path[i]].source;source;source=source->more)
      bundle_count[source->iseq]++;

  LOOP (i,seq->nsource_seq) {/* FOR EACH SEQ OVER THRESHOLD, ASSIGN bundle_id*/
/*    printf("bundle %d:\t%s\t%d/%d %d",bundle_id,seq->source_seq[i].name,
	   bundle_count[i],seq->source_seq[i].length,seq->source_seq[i].weight);*/
    if (seq->source_seq[i].bundle_id<0 /* NOT YET BUNDLED */
	&& seq->source_seq[i].length*minimum_fraction <= bundle_count[i]) {
/*      printf("   +++++++++++++++++");*/
      seq->source_seq[i].bundle_id = bundle_id; /* ASSIGN TO THIS BUNDLE */
      seq->source_seq[i].weight = 0; /* REMOVE FROM FUTURE heaviest_bundle */
      nseq_in_bundle++;
    }
  /*  printf("\n");*/
  }

  FREE(bundle_count);
  return nseq_in_bundle; /* RETURN COUNT OF SEQUENCES IN BUNDLE */
}




/** assigns weights for bundling based upon /hb_weight arguments
    in source_seq titles */

void assign_hb_weights(int nsource_seq,LPOSourceInfo_T source_seq[])
{
  int i,weight;
  char *p;
  LOOP (i,nsource_seq) {
    if (source_seq[i].title &&
	(p=strstr(source_seq[i].title,"/hb_weight="))) {
      weight=atoi(p+11);
      if (weight!=0){ /* 0 COULD MEAN atoi FAILED TO PARSE ARG.  IGNORE IT*/
	source_seq[i].weight = weight;
	fprintf(stderr,"assigned weight=%d to %s\n",source_seq[i].weight,source_seq[i].name);
      }
      else
	WARN_MSG(USERR,(ERRTXT,"hb_weight zero or unreadable: %s\nIgnored",p),"$Revision: 1.2 $");
    }
  }
}



/** generates the complete set of heaviest_bundle traversals of the the LPO
 seq, using iterative heaviest_bundle() and requiring that at least
 minimum_fraction of the positions in a sequence match the heaviest
 bundle path, for that sequence to be assigned to that bundle 
---------------------------------------------------------------
------------------------------------------------------------*/
void generate_lpo_bundles(LPOSequence_T *seq,float minimum_fraction)
{
  int nbundled=0,ibundle=0,path_length,iseq,count;
  LPOLetterRef_T *path=NULL;
  char name[256],title[1024];

  /*  assign_hb_weights(seq->nsource_seq,seq->source_seq); TURN THIS ON!!*/
  while (nbundled < seq->nsource_seq) {/* PULL OUT BUNDLES ONE BY ONE */
    path=heaviest_bundle(seq->length,seq->letter,/*GET NEXT HEAVIEST BUNDLE*/
			 seq->nsource_seq,seq->source_seq,&path_length);
    if (!path || path_length<10) /* ??!? FAILED TO FIND A BUNDLE ??? */
      goto premature_warning;
    sprintf(name,"CONSENS%d",ibundle);
    /* NEXT, MARK SEQUENCES THAT FIT THIS BUNDLE ADEQUATELY */
    count=assign_sequence_bundle_id(path_length,path,seq,ibundle,
				    minimum_fraction);
    sprintf(title,"consensus produced by heaviest_bundle, containing %d seqs",
	    count); /* DON'T INCLUDE CONSENSUS ITSELF IN THE COUNT! */
    iseq=add_path_sequence(path_length,path,seq,name,title);/*BUILD CONSENSUS*/
    seq->source_seq[iseq].bundle_id=ibundle++; /* INCREMENT BUNDLE ID */
    nbundled+=count; /* KEEP TRACK OF TOTAL SEQUENCES BUNDLED */

    if (count<1) {
    premature_warning:
      fprintf(stderr,"*** WARNING: bundling ended prematurely after %d bundles.\nNo sequences fit inside this last bundle.\nA total of %d sequences incuding consensus were bundled.\n\n",ibundle,nbundled);
      break;
    }
  }
}


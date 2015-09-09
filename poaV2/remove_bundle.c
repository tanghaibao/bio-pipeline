


#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"




int compact_links(LPOLetterLink_T *list,int old_to_new[])
{
  LPOLetterLink_T *link=NULL,*link_last=NULL,*next_link,*link_head=NULL;

  CALLOC(link,1,LPOLetterLink_T);
  memcpy(link,list,sizeof(LPOLetterLink_T));

  for (;link && link->ipos>=0;link=next_link){
    next_link=link->more;
    if (old_to_new[link->ipos]<0) /* THIS POSITION NO LONGER EXISTS! */
      free(link); /* DELETE THIS LINK ENTRY */
    else { /* COPY THIS BACK TO PREVIOUS LINK ENTRY: COMPACT THE LIST*/
      link->ipos = old_to_new[link->ipos]; /* REMAP TO NEW INDEX SYSTEM */
      if (link_last) /* CONNECT TO PREVIOUS NODE IN LIST */
	link_last->more=link;
      else /* THIS IS THE NEW HEAD OF THE LIST */
	link_head=link;
      link_last=link;
    }
  }
  if (link) /* AN EMPTY LINK I.E. link->ipos<0 ... JUNK IT */
    free(link);
  if (link_last) /* TERMINATE LAST NODE IN LIST */
    link_last->more=NULL;
  if (link_head) {
    memcpy(list,link_head,sizeof(LPOLetterLink_T));
    free(link_head);
    return 1; /* COMPACTED LINK LIST IS NON-EMPTY */
  }
  else { /* NOTHING LEFT IN LIST, SO BLANK IT */
    list->more=NULL;
    list->ipos= INVALID_LETTER_POSITION;
    return 0; /* COMPACTED LINK LIST IS EMPTY */
  }
}







int compact_sources(LPOLetterSource_T *list,int ibundle_delete,
		    LPOSourceInfo_T source_seq[])
{
  LPOLetterSource_T *source=NULL,*source_last=NULL,*next_source,*source_head=NULL;

  CALLOC(source,1,LPOLetterSource_T);
  memcpy(source,list,sizeof(LPOLetterSource_T));

  for (;source;source=next_source){
    next_source=source->more;
    if (source_seq[source->iseq].bundle_id == ibundle_delete)
      free(source); /* DELETE THIS SOURCE ENTRY */
    else { /* COPY THIS BACK TO PREVIOUS SOURCE ENTRY: COMPACT THE LIST*/
      if (source_last) /* CONNECT TO PREVIOUS NODE IN LIST */
	source_last->more=source;
      else /* THIS IS THE NEW HEAD OF THE LIST */
	source_head=source;
      source_last=source;
    }
  }
  if (source_last) /* TERMINATE LAST NODE IN LIST */
    source_last->more=NULL;
  if (source_head) {
    memcpy(list,source_head,sizeof(LPOLetterSource_T));
    free(source_head);
    return 1; /* COMPACTED SOURCE LIST IS NON-EMPTY */
  }
  else { /* NOTHING LEFT IN LIST, SO BLANK IT */
    list->more=NULL;
    list->ipos= INVALID_LETTER_POSITION;
    return 0; /* COMPACTED SOURCE LIST IS EMPTY */
  }
}






void reindex_compact_rings(LPOSequence_T *seq)
{
  int i,iring= -1,ring_start= -1;
  LOOPF (i,seq->length) { /* REMOVE ALL LINKS TO OLD, DELETED POSITIONS */
    if (iring==seq->letter[i].ring_id) {
      seq->letter[i].ring_id=ring_start; /* USE FIRST INDEX ON RING */
      seq->letter[i].align_ring= i-1; /* LINK TO LETTER TO LEFT */
    }
    else { /* START OF A NEW RING */
      if (ring_start>=0)
	seq->letter[ring_start].align_ring = i-1;
      iring=seq->letter[i].ring_id;
      seq->letter[i].ring_id=seq->letter[i].align_ring=ring_start=i;
    }
  }
  if (ring_start>=0)
    seq->letter[ring_start].align_ring = i-1;
}




#define DELETE_THIS_BUNDLE (-2)

int remove_bundle(LPOSequence_T *seq,int ibundle,int delete_all_others)
{
  int i,j=0,*old_to_new=NULL,new_length;

  if (delete_all_others) { /* INSTEAD OF DELETING THIS BUNDLE, */
    LOOP (i,seq->nsource_seq) /*MARK ALL OTHERS TO BE DELETED! */
      if (seq->source_seq[i].bundle_id!=ibundle)
	seq->source_seq[i].bundle_id = DELETE_THIS_BUNDLE;
    ibundle=DELETE_THIS_BUNDLE;
  }

  CALLOC(old_to_new,seq->length,int); /* CREATE MAPPING ARRAY */
  LOOPF (i,seq->length) {
    if (compact_sources(&seq->letter[i].source,ibundle,seq->source_seq)){
      if (i>j) /* COPY LETTER TO COMPACTED POSITION */
	memcpy(seq->letter+j,seq->letter+i,sizeof(LPOLetter_T));
      old_to_new[i]=j++; /* SAVE MAPPING FROM OLD TO NEW, COMPACTED POSITION*/
    }
    else {
      free_lpo_letters(1,seq->letter+i,FALSE); /* DUMP DATA FOR THIS LETTER */
      old_to_new[i]= INVALID_LETTER_POSITION;
    }
  }
  new_length=j;
  if (new_length<seq->length) /* ERASE UNUSED PORTIONS AFTER COMPACTED ARRAY*/
    memset(seq->letter+new_length,0,(seq->length - new_length)*sizeof(LPOLetter_T));

  LOOP (i,new_length) { /* REMOVE ALL LINKS TO OLD, DELETED POSITIONS */
    compact_links(&seq->letter[i].left,old_to_new);
    compact_links(&seq->letter[i].right,old_to_new);
  }
  seq->length=new_length;
  reindex_compact_rings(seq);

  FREE(old_to_new);
  return new_length;
}






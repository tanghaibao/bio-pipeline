

#include "default.h"
#include "seq_util.h"




void save_sequence_fields(Sequence_T *seq,
			  char seq_name[],char seq_title[],int length)
{
  STRNCPY(seq->name,seq_name,SEQUENCE_NAME_MAX);
  if (seq_title)
    seq->title=strdup(seq_title);
  else
    seq->title=strdup("untitled");
  seq->length=length; /* SAVE LENGTH */
}



int create_seq(int nseq,Sequence_T **p_seq,
	       char seq_name[],char seq_title[],char tmp_seq[],
	       int do_switch_case)
{
  int i,j;
  Sequence_T *seq;

  REBUFF(*p_seq,nseq,SEQUENCE_BUFFER_CHUNK,Sequence_T); /* ALLOCATE MEMORY*/
  seq= (*p_seq)+nseq; /* SET POINTER TO NEWLY ALLOCATED ELEMENT */

  for (i=j=0;tmp_seq[i];i++) /* ELIMINATE WHITE SPACE */
    if (!isspace(tmp_seq[i]))
      tmp_seq[j++]=tmp_seq[i];
  tmp_seq[j]='\0'; /* TERMINATE COMPRESSED STRING*/
  seq->sequence=strdup(tmp_seq); /* SAVE A DYNAMIC COPY */
  save_sequence_fields(seq,seq_name,seq_title,j);

  switch (do_switch_case) {
  case switch_case_to_lower:
    LOOP (i,seq->length)
      seq->sequence[i]=tolower(tmp_seq[i]);
    break;
  case switch_case_to_upper:
    LOOP (i,seq->length)
      seq->sequence[i]=toupper(tmp_seq[i]);
    break;
  }

  return 1;
}



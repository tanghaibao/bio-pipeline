
#include "default.h"


/*******************************************************************
  *
  *     STRINGPTR_CAT:
  *     holds string and count of last alloc'd size.  Conveniently
  *     manages static and dynamic strings, allowing auto-resize
  *     without worrying about memory management.
  *
  *     This function adds any string to the string ptr.  It even
  *     insulates against possibility that S2 is a substring of S1
  *     by pre-copying it.
  *
  ******************************************************************/

char *stringptr_cat_pos(stringptr *s1,const char s2[],int *pos) /* ~~g --- */
{
  int s2_len,total_len;
  char *stringptr_cat_temp_p=NULL,*s2_temp=NULL;
  
  if (s2 == NULL)					/* it might be better to make this a debug cope */
  	return s1->p;
  
  if (0 == s1->last_alloc) /* SAVE STATIC STRING */
    stringptr_cat_temp_p = s1->p; /* SAVE OLD */

  total_len=s2_len=strlen(s2)+1;
  CALLOC(s2_temp,s2_len,char); /* ALLOCATE TEMP STORAGE */
  memcpy(s2_temp,s2,s2_len); /* COPY THE STRING */

  if (s1->p) { /* CALCULATE ADDITIONAL SPACE NEEDED FOR ORIGINAL STRING */
    if (pos) /* USE THE CALLER-SUPPLIED STRING LENGTH */
      total_len += *pos;
    else /* OTHERWISE MEASURE IT */
      total_len += strlen(s1->p);
  }
  
  GETMEM(s1->p,total_len,
	 s1->last_alloc,STRINGPTR_BUFFER_CHUNK,char);

  if (stringptr_cat_temp_p) /* PUT OLD STATIC STRING BACK IN */ 
    strcpy(s1->p,stringptr_cat_temp_p); 

  if (pos) { /* ATTACH NEW STRING AT CALLER-SUPPLIED END-POSITION */
    strcpy(s1->p + *pos,s2_temp);
    *pos = total_len-1; /* EXCLUDE TERMINATOR */
  }
  else /* OTHERWISE strcat AS USUAL */
    strcat(s1->p,s2_temp);

  FREE(s2_temp);  /* FREE THE TEMP STORAGE */

  return s1->p;
}


char *stringptr_cat(stringptr *s1,const char s2[]) /* ~~g --- */
{
  return stringptr_cat_pos(s1,s2,NULL);
}


char *stringptr_cpy(stringptr *s1,const char s2[]) /* ~~g --- */
{
  GETMEM(s1->p,strlen(s2)+1,s1->last_alloc,STRINGPTR_BUFFER_CHUNK,char);
  strcpy(s1->p,s2);

  return s1->p;
}






/**stringptr_free*************************************************
  *
  *     stringptr_free:
  *     AUTHOR: tal
  *     Wed Jul 27 03:04:41 PDT 1994
  *
  *     Frees stringptr type.
  *
  ***************************************************************/

int stringptr_free(stringptr *s) /* ~~g --- */
{
  FREE(s->p);
  s->last_alloc=0;
  return 0;
}



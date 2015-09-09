#ifndef DEFAULT_HEADER_INCLUDED
#define DEFAULT_HEADER_INCLUDED 1

#ifndef MODULE_NAME
#define MODULE_NAME "main"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#include "black_flag.h"



typedef void *voidptr;  /* ~~e: should be moved out to generic typing header
			   --- */
typedef int (*funptr)();

#define LOOPB(i,size) for ((i)=(size);(i)-- >0;)
#define LOOP(i,size) for ((i)=(size);(i)-- >0;)
#define LOOPF(i,size) for ((i)=0;(i)<(size);(i)++)
#define LOOP_FINISHED(i,size) ((i)<0 || (i)>=(size))


/**@memo example usage of argument reading macros 
  FOR_ARGS(i,argc) {
    ARGMATCH_VAL("-tolower",do_switch_case,switch_case_to_lower);
    ARGMATCH("-seq_err_report",count_sequence_errors);
    ARGGET("-printmatrix",print_matrix_letters);
    NEXTARG(matrix_filename);
  } */

#define FOR_ARGS(INDEX,ARGC) for (INDEX=1;INDEX<ARGC;INDEX++)
#define ARGMATCH(FLAGSTR,VAR) if (0==strcmp(argv[i],FLAGSTR)) (VAR)=1
#define ARGMATCH_VAL(FLAGSTR,VAR,VAL) if (0==strcmp(argv[i],FLAGSTR)) (VAR)=(VAL)
#define ARGGET(FLAGSTR,VAR) if (0==strcmp(argv[i],FLAGSTR)) {(VAR)=argv[++i];continue;}
#define NEXTARG(VAR) if ('-'!=argv[i][0] && !(VAR)) {(VAR)=argv[i];continue;}




#ifndef TRUE  /* DEFINE true AND false SYMBOLS IN CASE NOT ALREADY THERE*/
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif
  
#if defined(__STDC__) || defined(__ANSI_CPP__)  /*###################*/
/* ~~a: ansi-dependent concatenation macro */
#define CONCAT_MACRO0(a,b)  a ## b
#else
#define CONCAT_MACRO0(a,b) a/**/b
#endif /* !__STDC__ ##################################################*/
  
#define CONCAT_MACRO(a,b)  CONCAT_MACRO0(a,b)

#if defined(__STDC__) || defined(__ANSI_CPP__)  /*###################*/
#define STRINGIFY(NAME) # NAME
#else
#define STRINGIFY(NAME) "NAME"
#endif /* !__STDC__ ##################################################*/

/* #define PRINT_DEBUG(MESSAGE)  printf MESSAGE */


/***********************************************************************
  *
  *     STRNCPY: guaranteed overflow-protected, null-terminated strcpy
  *
  *     NB: this macro forces termination at end position S1[LEN-1]
  *     which will cause non-intuitive side-effects ESPECIALLY
  *     if LEN is incorrect!!!!  ie. STRNCPY will garbage your
  *     memory by sticking zeroes into it if you give it the wrong
  *     LEN.
  *
  *     Make sure that the LEN you give it is correct!!!!!!!
  *
  **********************************************************************/

#define STRNCPY(S1,S2,LEN) \
  (strncpy(S1,S2,LEN),(((char *)(S1))[(LEN)-1]='\0'))








/*******************************************************************
  *
  *     STRINGPTR:
  *     holds string and count of last alloc'd size.  Conveniently
  *     manages static and dynamic strings, allowing auto-resize
  *     without worrying about memory management.
  *
  ***********************************************************************/

typedef struct { /* ~~d */
  char *p;
  int last_alloc;
} stringptr; /* --- */


#define STRINGPTR_BUFFER_CHUNK 4096
#define STRINGPTR_EMPTY_INIT {NULL,0}

char *stringptr_cat(stringptr *s1,const char s2[]);
char *stringptr_cpy(stringptr *s1,const char s2[]);
int stringptr_free(stringptr *s);






#ifdef DEBUG
#define PRINT_DEBUG(CODE,MESSAGE) {if (CODE<=2) fprintf MESSAGE ;}
#else
#define PRINT_DEBUG(CODE,MESSAGE) {if (CODE<=1) fprintf MESSAGE ;}
#endif




#define MEM0_REQUEST_BAD 0

/* DEFAULT ACTION TO PERFORM IF malloc FAILS */
#define MALLOC_FAILURE_ACTION abort()
/* DEFAULT ACTION TO PERFORM IF realloc FAILS */
#define REALLOC_FAILURE_ACTION abort()



#define CALLOC(memptr,N,ATYPE) \
/*  printf("CALLOC(%s=%d,%s=%d,sizeof(%s)=%d), %s %d\n",STRINGIFY(memptr),(memptr), \
	 STRINGIFY(N),(N),STRINGIFY(ATYPE),sizeof(ATYPE),__FILE__,__LINE__);*/ \
  if ((N)<=0)   {                     \
    if (MEM0_REQUEST_BAD)  {\
      fprintf(stderr,"%s, line %d: *** invalid memory request: %s[%d].\n",\
	      __FILE__,__LINE__,STRINGIFY(memptr),(N));   \
      MALLOC_FAILURE_ACTION;                                            \
    }\
  }                                                                 \
  else if (NULL == ((memptr)=(ATYPE *)calloc((size_t)(N),sizeof(ATYPE))))  { \
    fprintf(stderr,"%s, line %d: *** out of memory \n",__FILE__,__LINE__);                \
    fprintf(stderr,"Unable to meet request: %s[%d]\n",STRINGIFY(memptr),(N));    \
    fprintf(stderr,"requested %d x %d bytes \n",(N),sizeof(ATYPE));   \
    MALLOC_FAILURE_ACTION;                                            \
  }









/*******************************************************************
  *
  *     FREE:
  *     dumps storage if a non-zero pointer, and resets pointer to
  *     NULL to indicate its data freed.
  *
  ******************************************************************/

#define FREE(P) if (P) {free(P);(P)=NULL;}


#define REALLOC(memptr,NUM,ATYPE) \
/*  printf("REALLOC(%s=%d,%s=%d,sizeof(%s)=%d), %s %d\n",STRINGIFY(memptr), \
    (memptr),STRINGIFY(NUM),(NUM),STRINGIFY(ATYPE),sizeof(ATYPE),\
    __FILE__,__LINE__);*/ \
  if ((NUM)<=0)   {                            \
    if (MEM0_REQUEST_BAD)  {\
      fprintf(stderr,"%s, line %d: *** invalid memory request: %s[%d].\n",\
	      __FILE__,__LINE__,STRINGIFY(memptr),(NUM));   \
      REALLOC_FAILURE_ACTION;                                             \
    }\
  }                                                                 \
  else { \
    voidptr temp_ptr=realloc((void *)(memptr),sizeof(ATYPE)*(size_t)(NUM));\
    if (temp_ptr) \
      (memptr)=(ATYPE *)temp_ptr;\
    else { \
      fprintf(stderr,"%s, line %d: *** out of memory \n",__FILE__,__LINE__); \
      fprintf(stderr,"Unable to meet request: %s\n",STRINGIFY(memptr));  \
      fprintf(stderr,"requested %d x %d bytes \n",(NUM),sizeof(ATYPE));   \
      REALLOC_FAILURE_ACTION;                                             \
    } \
  }


/**********************************************************************
  *
  *     REBUFF
  *     GUARANTEES ALLOCATION OF NULL-INITIALIZED BLOCK OF BUF ELEMENTS
  *     OF THE CORRECT DATA TYPE & SIZE, SUFFICIENT TO HOLD AT LEAST ONE
  *     ELEMENT ADDED TO THE CURRENT SIZE OF THE ARRAY
  *
  **********************************************************************/

#define REBUFF(memptr,NUM,BUF,ATYPE) if ((NUM)<=0) { \
    CALLOC(memptr,BUF,ATYPE) } \
  else if (0 == (NUM)%(BUF)) { \
    char *p_rebuff_size;\
    REALLOC(memptr,(NUM)+(BUF),ATYPE); \
    p_rebuff_size=(char *)(memptr)+(NUM)*sizeof(ATYPE); \
    memset((voidptr)p_rebuff_size,0,(size_t)((BUF)*sizeof(ATYPE)));\
  }









/*******************************************************************
  *
  *     GETMEM:
  *     allocates memory either via malloc or realloc, guaranteeing
  *     space to store at least one more data entry.  Data is
  *     alloc'd in blocks of size BUF; the total # alloc'd is
  *     kept in LAST.  The algorithm can go wrong if LAST
  *     is incorrect.
  *
  ******************************************************************/

#define GETMEM(memptr,NUM,LAST,BUF,TYPE) \
  if ((LAST)<=0) {\
    CALLOC(memptr,((NUM)+1)-((NUM)+1)%(BUF)+(BUF),TYPE);\
    (LAST)=((NUM)+1)-((NUM)+1)%(BUF)+(BUF);\
  }\
  else if (((NUM)+1)>(LAST))  { \
    REALLOC(memptr,((NUM)+1)-((NUM)+1)%(BUF)+(BUF),TYPE);\
    (LAST)=((NUM)+1)-((NUM)+1)%(BUF)+(BUF);\
  }





#endif

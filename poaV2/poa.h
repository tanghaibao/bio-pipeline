

#ifndef POA_HEADER_INCLUDED
#define POA_HEADER_INCLUDED




/** MAXIMUM GAP LENGTH TRACKED IN align_lpo;
    LARGER THAN THIS WILL BE CAPPED
    (DEFAULT VALUE) */
#ifndef TRUNCATE_GAP_LENGTH
#define TRUNCATE_GAP_LENGTH 16
#endif

/** LENGTH OVER WHICH GAP PENALTY DECAYS IN align_lpo
    (DEFAULT VALUE) */
#ifndef DECAY_GAP_LENGTH
#define DECAY_GAP_LENGTH 0
#endif

#define REV_COMP_STRING "/rev_comp"


/** THE NULL LETTER-REFERENCE */
#define INVALID_LETTER_POSITION (-1)

typedef int LPOLetterRef_T;


typedef int LPOScore_T;

 /** NEEDED FOR seq_util.h */
typedef LPOScore_T ResidueScore_T;
#define RESIDUE_SCORE_DEFINED

/** linked list for storing source origin (sequence position)
 from which this letter was derived */
struct LPOLetterSource_S {
 /** index of the sequence, referencing the source_seq[] array*/
  int iseq;
 /** index of the corresponding position in that sequence */
  LPOLetterRef_T ipos;
 /** next node in the linked list */
  struct LPOLetterSource_S *more;
} ;

typedef struct LPOLetterSource_S LPOLetterSource_T;


/** linked list for connecting an LPOLetter to either right
 or left */
struct LPOLetterLink_S {
 /** ADJACENT LETTER LINKED TO THIS LETTER */	
  LPOLetterRef_T ipos;
#ifdef USE_WEIGHTED_LINKS
 /** transition cost for traversing this link */
  LPOScore_T score;
#endif
 /** next node in the linked list */
  struct LPOLetterLink_S *more;
} ;

typedef struct LPOLetterLink_S LPOLetterLink_T;


/** the chunk size for allocating additional
    letters in an LPOLetter_T array */
#define LPO_LETTER_BUFFER_CHUNK 64


/** Structure for storing individual LPO Letters*/
struct LPOLetter_S {
 /** ADJACENT LETTER(S) TO THE LEFT */	
  LPOLetterLink_T left;
 /** ADJACENT LETTER(S) TO THE RIGHT */
  LPOLetterLink_T right;
 /** SOURCE SEQ POSITION(S) */
  LPOLetterSource_T source;
 /** CIRCULAR LIST OF ALIGNED POSITIONS */
  LPOLetterRef_T align_ring;
 /** MINIMUM INDEX OF ALL POSITIONS ON THE RING */
  LPOLetterRef_T ring_id;
  /** SCORE FOR BALANCING PARTIAL ORDER EFFECTS ON MATRIX NEUTRALITY */
  float score;
 /** THE ACTUAL RESIDUE CODE! */
  char letter;
} ;



typedef struct LPOLetter_S LPOLetter_T;


/** maximum length of a sequence name */
#define SEQUENCE_NAME_MAX 32
/** buffer chunk size for expanding a block of seq storage */
#define SEQUENCE_BUFFER_CHUNK 8

/** buffer chunk size for expanding a source_seq[] array */
#define SOURCE_SEQ_BUFFER_CHUNK 16



#define NUMDATA_BUFFER_CHUNK 4
/** storage for quantitative data attached to a sequence */
struct LPONumericData_S { /** */
  char name[SEQUENCE_NAME_MAX]; /** */
  char *title; /** */
  double *data;
};
typedef struct LPONumericData_S LPONumericData_T;


/** Structure for storing individual source sequence information,
 stuff like name, title etc. */
struct LPOSourceInfo_S { /** */
  char name[SEQUENCE_NAME_MAX]; /** */
  char *title; /** */
  char *sequence; /** */
  int *seq_to_po; /** */
  int *po_to_seq; /** */
  LPONumericData_T *data; /** */
  int ndata; /** */
  int length; /** */
  int istart;
 /** FOR PURPOSES OF HEAVIEST BUNDLE CALCULATION */
  int weight;
 /** WHAT BUNDLE IS THIS A MEMBER OF? */
  int bundle_id;
};

typedef struct LPOSourceInfo_S LPOSourceInfo_T;

/** the NULL bundle-reference */
#define NO_BUNDLE (-1)

/** bundle-reference meaning "include all bundles" */
#define ALL_BUNDLES (-1)


/** holder for an LPO sequence, its letters, 
  and associated information */
struct LPOSequence_S {/** */
  int length;/** */
  LPOLetter_T *letter;/** */
  char *title;/** */
  char *sequence;/** */
  char name[SEQUENCE_NAME_MAX];/** */
  int nsource_seq;/** */
  LPOSourceInfo_T *source_seq;
};

typedef struct LPOSequence_S LPOSequence_T;


typedef LPOSequence_T Sequence_T;


/**@memo GENERAL FORM IS seq_y[j].left.ipos */
#define SEQ_Y_LEFT(j) (j-1)
#define SEQ_Y_RIGHT(j) (j+1)



/**@memo Data structure for analyzing sequence differences in MSA*/
struct LPOLetterCount_S {
  unsigned int is_error:2;
  unsigned int meets_criteria:1;
  unsigned int seq_count:29;
};

typedef struct LPOLetterCount_S LPOLetterCount_T;

/** classification of sequence differences */
enum {
  no_error,
  substitution_error,
  insertion_error,
  deletion_error,
  max_error_states
};







/** DON'T ALLOCATE MORE THAN THIS TOTAL AMOUNT OF MEMORY 
---------------------------------------------------------------
---------------------------------------------------------------
*/
#define POA_MAX_ALLOC 300000000


#endif






/**@name The lpo library*/
/*@{*/
/**@memo This set of web pages documents the functionality of the lpo 
function library.  This is a set of C functions for reading, writing,
creating, manipulating, and aligning partial order sequences.  These 
functions divide into several groups:
\begin{itemize} 
\item \URL[File utilities]{General.html#read_fasta}: reading and writing 
FASTA and po files
\item \URL[lpo utilities]{General.html#add_lpo_link}: creating, fusing, 
freeing, manipulating lpo data
\item \URL[alignment]{General.html#align_lpo}: aligning one or more linear 
sequences to an lpo
\item analysis: analyzing lpo structure, e.g. to find consensus
\end{itemize}
*/

/**@memo Click \URL[here]{../poa} for more information about partial 
order alignment.*/

/*@}*/



/**@name linking to the lpo library */
/*@{*/
/**@memo To use function from this library in your code, you must 
 do two things.  First you must 
include <lpo.h> in your source files, to access the prototypes.  
 Second, when you compile, you must tell the compiler where the 
 lpo header and library files are located.  {\bfNB it appears gcc 
loads libraries in reverse order of the command line arguments, so 
you have to specify your source files BEFORE the -llpo library argument 
on the command line, or the linker will give you unresolved reference 
errors}.  e.g.  \begin{verbatim}

gcc -o myprog myfile.c -I~leec/lib/include -L~leec/lib -llpo
\end{verbatim}
*/

/**@memo Click \URL[here]{../poa} for more information about partial 
order alignment.*/

/*@}*/






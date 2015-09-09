
#include "default.h"
#include "seq_util.h"



/** reads FASTA formatted sequence file, and saves the sequences to
  the array seq[]; any comment line preceded by a hash-mark will be saved
  to comment */
int read_fasta(FILE *seq_file,Sequence_T **seq,
	       int do_switch_case,char **comment)
{
  int c,nseq=0,length=0;
  char seq_name[FASTA_NAME_MAX]="",
  line[SEQ_LENGTH_MAX],seq_title[FASTA_NAME_MAX]="";
  char *p;
  stringptr tmp_seq=STRINGPTR_EMPTY_INIT;
	
 /* read in sequences */
  while (fgets(line,sizeof(line)-1,seq_file)) {
    if ((p=strrchr(line,'\n'))) /* REMOVE NEWLINE FROM END OF LINE */
      *p= '\0'; /* TRUNCATE THE STRING */
    switch (line[0]) {
    case '#':  /* SEQUENCE COMMENT, SAVE IT */
      if (comment) /* SAVE COMMENT FOR CALLER TO USE */
	*comment = strdup(line+1);
      break;

    case '>':  /* SEQUENCE HEADER LINE */
      if (seq_name[0] && tmp_seq.p && tmp_seq.p[0]) { /* WE HAVE A SEQUENCE, SO SAVE IT! */
	if (create_seq(nseq,seq,seq_name,seq_title,tmp_seq.p,do_switch_case))
	  nseq++;
      }
      seq_name[0]='\0';
      if (sscanf(line+1,"%s %[^\n]",  /* SKIP PAST > TO READ SEQ NAME*/
	     seq_name,seq_title)<2)
	strcpy(seq_title,"untitled"); /* PROTECT AGAINST MISSING NAME */
      if (tmp_seq.p)
	tmp_seq.p[0]='\0'; /* RESET TO EMPTY SEQUENCE */
      length=0;
      break;

    case '*': /* IGNORE LINES STARTING WITH *... DON'T TREAT AS SEQUENCE! */
      break;

    default:  /* READ AS ACTUAL SEQUENCE DATA, ADD TO OUR SEQUENCE */
      if (seq_name[0]) /* IF WE'RE CURRENTLY READING A SEQUENCE, SAVE IT */
	stringptr_cat_pos(&tmp_seq,line,&length);
    }

    c=getc(seq_file); /* ?FIRST CHARACTER IS UNIGENE CLUSTER TERMINATOR? */
    if (c==EOF)
      break;
    else {
      ungetc(c,seq_file); /* PUT THE CHARACTER BACK */
      if (c=='#' && nseq>0) /* UNIGENE CLUSTER TERMINATOR, SO DONE!*/
	break;
    }
  }
  if (seq_name[0] && tmp_seq.p && tmp_seq.p[0]) { /* WE HAVE A SEQUENCE, SO SAVE IT! */
    if (create_seq(nseq,seq,seq_name,seq_title,tmp_seq.p,do_switch_case))
      nseq++;
  }
  stringptr_free(&tmp_seq);
  return nseq; /* TOTAL NUMBER OF SEQUENCES CREATED */
}

/**@memo example: reading FASTA format file: 
    seq_ifile=fopen(seq_filename,"r");
    if (seq_ifile) {
      nseq=read_fasta(seq_ifile,&seq,do_switch_case,&comment);
      fclose(seq_ifile);
    }
*/



/** writes a FASTA formatted file, saving the sequence given in seq[] */
void write_fasta(FILE *ifile,char name[],char title[],char seq[])
{
  int j;

  fprintf(ifile,">%s %s\n",name,title? title : "untitled");
  for (j=0;j<strlen(seq);j+=60)
    fprintf(ifile,"%.60s\n",seq+j);

  return;
}


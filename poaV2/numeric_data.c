

#include <lpo.h>

/**@memo create a new LPONumericData_T record for the designated source_seq.  
 Dynamically allocates data storage array equal in length to the length of 
 the source_seq.  Initializes the array values to initial_value if non-zero.  
 Returns a pointer to the new LPONumericData_T and also adds it to the 
 source_seq->data[] list.
*/
LPONumericData_T *new_numeric_data(LPOSourceInfo_T *source_seq,
				   char name[],
				   char title[],
				   double initial_value)
{
  int i;
  LPONumericData_T *data=NULL;
  REBUFF(source_seq->data,source_seq->ndata,NUMDATA_BUFFER_CHUNK,LPONumericData_T);
  data=source_seq->data + source_seq->ndata++;

  STRNCPY(data->name,name,SEQUENCE_NAME_MAX); /* COPY NAME AND TITLE */
  if (title)
    data->title=strdup(title);

  CALLOC(data->data,source_seq->length,double); /* ALLOCATE THE ARRAY */
  if (initial_value) /* INITIALIZE VALUES IF DESIRED */
    LOOP (i,source_seq->length)
      data->data[i]=initial_value;
  return data; /* RETURN POINTER TO THE NEW DATA HOLDER */
}



LPONumericData_T *cp_numeric_data(LPOSourceInfo_T *source_seq,
				  LPONumericData_T *data)
{
  int i;
  LPONumericData_T *new_data;
  new_data=new_numeric_data(source_seq,data->name,data->title,0);
  LOOP (i,source_seq->length) /* COPY ALL THE VALUES */
    new_data->data[i]=data->data[i];
  return new_data;
}



/**@memo finds LPONumericData from the source_seq, matching the specified 
   name, or returns NULL if not found. */
LPONumericData_T *find_numeric_data(LPOSourceInfo_T *source_seq,
				   char name[])
{
  int i;
  LOOP (i,source_seq->ndata)
    if (0==strcmp(source_seq->data[i].name,name)) /*FOUND IT. RETURN POINTER*/
      return source_seq->data+i;
  return NULL; /* NOT FOUND! */
}



/**@memo frees the set of numeric_data passed as arguments.  If requested, 
 will also free the block of memory for the array of entries data[]. */
void free_lpo_numeric_data(int ndata,LPONumericData_T *data,
			   int please_free_block)
{
  int i;
  LOOP (i,ndata) { /* DUMP ASSOCIATED ARRAYS */
    FREE(data[i].title);
    FREE(data[i].data);
  }
  if (please_free_block) /* DUMP THE BLOCK ITSELF */
    free(data);
}





/**@memo creates one or more new numeric_data for a given sequence, 
  based on the presence of corresponding named numeric_data.  Specifically, 
  the list of set_names[] is processed one by one, creating a source_name
  according to the source_name_fmt string, and finding a numeric_data 
  entry with that name.  If it is not found, the routine calls exit(-1).  
  If it is found, a new set of numeric_data is created, with a name 
  created according to target_name_fmt, and titled according to 
  the title_fmt string and source_data->title. */
void new_numeric_data_sets(LPOSourceInfo_T *source_seq,
			   int nset,char *set_names[],
			   char source_name_fmt[],
			   char target_name_fmt[],
			   char title_fmt[])
{
  int j;
  LPONumericData_T *data;
  char name[256],title[4096];

  LOOPF (j,nset) {
    sprintf(name,source_name_fmt,set_names[j]); /* GENERATE NAME TO MATCH*/
    data=find_numeric_data(source_seq,name); /*FIND SOURCE DATA*/
    if (!data) {
      WARN_MSG(USERR,(ERRTXT,"*** could not find dataset %s for seq %s.\nExiting\n\n",
	      name,source_seq->name),"$Revision: 1.2 $");
      exit(-1);
    }
    sprintf(name,target_name_fmt,set_names[j]); /* GENERATE NEW NAME */
    sprintf(title,title_fmt,data->title);
    new_numeric_data(source_seq,name,title,0.); /* CREATE NEW ARRAY */
  }
}




/**@memo reads a stream of FASTA-formatted numeric data, and stores them in 
the corresponding set of source_seq, matching the sequences by name.  
Multiple numeric data entries can be read from a single stream. */
void read_numeric_data(int nsource_seq,
		       LPOSourceInfo_T source_seq[],
		       FILE *ifile)
{
  int i,j;
  char line[4096],seq_name[128],data_name[1024],title[2048];
  LPONumericData_T *data;

  while (fgets(line,sizeof(line),ifile)) {
    title[0]='\0';
    if (sscanf(line,">%s NUMERIC_DATA=%s %s",seq_name,data_name,title)>=2) {
      LOOP (i,nsource_seq) /* FIND THE MATCHING SEQUENCE */
	if (0==strcmp(seq_name,source_seq[i].name)) /* MATCH */
	  break;
      if (LOOP_FINISHED(i,nsource_seq)) /* SEQ NOT FOUND!! */
	WARN_MSG(USERR,(ERRTXT,"Error! NUMERIC_DATA %s, sequence %s does not exist.  Skipping.\n\n",data_name,seq_name),"$Revision: 1.2 $");
      else { /* FOUND THE SEQ, SAVE THE DATA */
	if (data=find_numeric_data(source_seq+i,data_name)) /*REUSE EXISTING*/
	  WARN_MSG(WARN,(ERRTXT,"NUMERIC_DATA %s already exists on sequence %s.  Overwriting.\n",data_name,seq_name),"$Revision: 1.2 $");
	else /* CREATE A NEW DATA HOLDER */
	  data=new_numeric_data(source_seq+i,data_name,title,0.);
	LOOPF (j,source_seq[i].length) /* READ IN THE VALUES */
	  fscanf(ifile," %lf",data->data+j);
      }
    }
  }
}



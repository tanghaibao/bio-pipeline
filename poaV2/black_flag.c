


#include "default.h" /* ~~I */




char ERRTXT[1024]
  ="";
char *Program_name="black_flag";
char *Program_version="unknown";
int Already_reported_crash=0;

int black_flag(int bug_level,
	       char sourcefile[],
	       int sourceline,
	       char sourcefile_revision[])
{
  static int last_line= -1;
  char *error_names[max_black_flag_type]
    ={"CRASH","DIED","EXCEPTION","BAD_DATA","WARNING","DEBUG"};

  switch (bug_level) {
#ifndef DEBUG_USER_VERSION
  case DEBUG_black_flag_type:  /* IN DEV'T VERSION JUST CRASH!!! */
#ifdef DEBUG_VERSION
  case TRAP_black_flag_type: /* FOR DEBUG VERSION AND TRAP, CAUSE A CORE DUMP*/
#endif
#endif
  case CRASH_black_flag_type:  /* COULD INCLUDE MECHANISMS TO SEND EMAIL? */
    /* print message at level 1, i.e. if we are printing anything at all */
    PRINT_DEBUG(1,(DBOUT,"black_flag: %s %s:%s %s %s,%d\n%s\nend_black_flag\n",
		   error_names[bug_level],
		   Program_name,Program_version,
		   sourcefile,sourcefile_revision,sourceline,
		   ERRTXT[0]? ERRTXT:""));
    if (Already_reported_crash)
      return 1;
    Already_reported_crash=1;
    abort();  /* FORCE A CORE DUMP */


  default: /* JUST PRINT THE ERROR */
    PRINT_DEBUG(1,(DBOUT,"black_flag: %s %s:%s %s %s,%d\n%s\nend_black_flag\n",
		   error_names[bug_level],
		   Program_name,Program_version,
		   sourcefile,sourcefile_revision,sourceline,
		   ERRTXT[0]? ERRTXT:""));
    break;
    
  }

  ERRTXT[0]='\0'; /* RESET THE ERROR TEXT */
  last_line=sourceline;
  return 1;  /* SEND SIGNAL TO HANDLER CLAUSE TO DEAL WITH THIS ERROR */
}








void handle_crash(int sigcode)
{
  char *crash_mode;
  
  if (Already_reported_crash)
    exit(-1); /* IN ENDLESS LOOP REPORTING CRASH OVER & OVER? */

  Already_reported_crash=1;
  black_flag(CRASH_black_flag_type,"",0,"");
  if ((crash_mode=getenv("HANDLE_CRASH")) && 0==strcmp(crash_mode,"NOCORE"))
    exit(-1);
  else
    return; /*ENVIRONMENT SETTING ASKS US TO DUMP A CORE IMMEDIATELY */
}





int handle_crash_init(void (*crash_fun)())
{
#define HANDLE_CRASH_MAX 5
  int i,signal_type[HANDLE_CRASH_MAX]
    ={SIGSEGV,
#ifdef SIGBUS /* LINUX DOESN'T HAVE BUS ERRORS? */
	SIGBUS,
#else
	SIGSEGV, /* REUSE SEGV AS DUMMY ENTRY */
#endif
	SIGABRT,SIGFPE,SIGTRAP}; /*LIST OF SIGNALS TO HANDLE*/
  if (!crash_fun) /* NO CRASH FUNCTION SUPPLIED, SO RESET TO DEFAULT */
    crash_fun = SIG_DFL; /* RESET TO STANDARD CRASH BEHAVIOR */

  LOOP (i,HANDLE_CRASH_MAX) /* SET THIS HANDLER FOR ALL OUR SIGNALS */
    signal(signal_type[i],crash_fun);
  return 0;
}

void black_flag_init_args(int narg,char *arg[],char progversion[])
{
  int i,len=0;
  LOOP (i,narg)
    len+=strlen(arg[i])+2;
  CALLOC(Program_name,len,char);
  LOOPF (i,narg) {
    strcat(Program_name,arg[i]);
    strcat(Program_name," ");
  }
  if (progversion)
    Program_version=progversion;
  handle_crash_init(handle_crash);
}



void black_flag_init(char progname[],char progversion[])
{
  black_flag_init_args(1,&progname,progversion);
}


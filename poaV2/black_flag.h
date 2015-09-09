


#ifndef BLACK_FLAG_HEADER_INCLUDED
#define BLACK_FLAG_HEADER_INCLUDED 1

#include <signal.h>

extern char ERRTXT[];
#define DBOUT stderr

enum {
  CRASH_black_flag_type,
  TRAP_black_flag_type,
  COPE_black_flag_type,
  USERR_black_flag_type,
  WARN_black_flag_type,
  DEBUG_black_flag_type,
  max_black_flag_type
  };


#define NOERRMSG (ERRTXT,"")





/********************************************************************
  *
  *     IF_PARANOID:
  *     IF_PARANOID(CONDITION,REVISION,MESSAGE)
  *
  *     error checking that will significantly slow execution or
  *     otherwise is desirable only in situations of extreme
  *     unction, e.g. during our debugging!!!
  *
  *     The beauty of the IF_PARANOID idea is that you should sprinkle
  *     it liberally EVERYWHERE in your code without thought for
  *     whether it is necessary or might hurt performance,
  *     because it will NOT even be compiled into the program in
  *     the production version!!!!!!
  *
  *     use IF_PARANOID checks EVERYWHERE you can think of
  *     definite error signals, even if you think "That shouldn't
  *     EVER happen".
  *
  *     if the IF_PARANOID CONDITIONAL is TRUE, the program will abort();
  *
  *
  *
  *     Also, our emacs custom highlighting system has been programmed
  *     to show IF_PARANOID and IF_DEBUG lines in a dim gray, so
  *     that these extensive error checks will not visually obscure
  *     the layout & organization of your algorithms.
  *
  *     black_flag() is smart about printing both version information
  *     such as the vmake version name the executable was created by,
  *     and also the exact revision number of the file in which the
  *     error occured.
  *
  *
  *-------------------------------------------------
  *
  *     IF_DEBUG:
  *     IF_DEBUG(CONDITION,REVISION,MESSAGE)
  *
  *     also not expected to be included in a final production
  *     release, but should not impact performance so significantly that
  *     it's unpleasant to test such a version in regular use patterns.
  *     Classic examples would be fairly pedantic checks at the entry
  *     and exit of all functions, testing for conditions "that shouldn't
  *     ever happen."
  *
  *     if the IF_DEBUG CONDITIONAL is TRUE, the program will abort();
  *
  *
  *-------------------------------------------------
  *
  *     IF_GUARD:
  *     IF_GUARD(CONDITION,REVISION,MESSAGE,LEVEL)
  *
  *     checks INCLUDED in final production versions, but taking
  *     advantage of the black_flag system to allow us to easily
  *     control what will be done in response to an error
  *     in any given executable, via compile-time flags:
  *
  *     e.g.
  *     print error info on stderr, or to special log files,
  *
  *     force a core dump,
  *
  *     run dbx via an auto script to generate a stack frame,
  *     and mail the results to develop-support@mag.com,
  *     etc.
  *
  *     IF_GUARD differs from IF_PARANOID and IF_DEBUG in that it
  *     requires an error_level argument, which must be one of
  *
  *     CRASH                     ... the error is fatal, abort()
  *
  *     TRAP                      ... the error is being trapped, but
  *                                   not handled. e.g. quiting from
  *                                   a function because some essential
  *                                   file was missing...
  *
  *     COPE                      ... the error is being handled nicely.
  *                                   the handler code following IF_GUARD
  *                                   knows how to correct for the situation.
  *
  *     USERR                     ... the user appears to have done something
  *                                   that makes no sense; they must be
  *                                   confused by our interface.
  *
  *     WARN                      ... suspicious data or input, not
  *                                   definitely an error.
  *
  *
  *
  *-------------------------------------------------
  *
  *     Examples:
  * 
  *     IF_PARANOID((iatom<0),4.6,(ERRTXT,"wacky iatom=%d",iatom));
  *
  *     IF_DEBUG((iatom<0),4.6,(ERRTXT,"wacky iatom=%d",iatom));
  *
  *     IF_GUARD((iatom<0),4.6,(ERRTXT,"wacky iatom=%d",iatom),CRASH);
  *
  *     IF_GUARD((iatom<0),4.6,(ERRTXT,"wacky iatom=%d",iatom),COPE) {
  *       put some code to handle the error condition here;
  *     }
  *
  *     USE THE 4.6 KEYWORD TO PUT IN VERSION NUMBERS AUTOMATICALLY.
  *
  *******************************************************************/


#if defined(PARANOID_VERSION)  /*???????????????????*/
#ifndef DEBUG_VERSION
#define DEBUG_VERSION
#endif

#define IF_PARANOID(CONDITION,REVISION,MESSAGE) \
  if (CONDITION) {\
    sprintf MESSAGE ;\
    black_flag(DEBUG_black_flag_type,__FILE__,__LINE__,STRINGIFY(REVISION));\
  }

#else /*????????????????????????????????????????????*/
#define IF_PARANOID(CONDITION,REVISION,MESSAGE)
#endif  /* !PARANOID_VERSION  ??????????????????????*/










#if defined(DEBUG_VERSION) /*??????????*/
#define IF_DEBUG(CONDITION,REVISION,MESSAGE) \
  if (CONDITION) {\
    sprintf MESSAGE ;\
    black_flag(DEBUG_black_flag_type,__FILE__,__LINE__,STRINGIFY(REVISION));\
  }

#else  /*???????????????????????????????????????????????????????????*/
#define IF_DEBUG(CONDITION,REVISION,MESSAGE)
#endif  /* !DEBUG_VERSION   ????????????????????????????????????????*/








/*********************************************************************
  *
  *     IF_GUARD:
  *     the basic black_flag macro, for production version trapping
  *     and handling of errors.
  *
  *     Essentially adds the flexibility of handling / recording
  *     errors however you like in black_flag().  Also, the programmer
  *     can provide, or omit, a clause following this macro that
  *     will only be executed if the CONDITION is true, allowing
  *     any kind of handling you wish.
  *
  ********************************************************************/


#define IF_GUARD(CONDITION,REVISION,MESSAGE,LEVEL) \
  if ((CONDITION) ? \
      (sprintf MESSAGE,\
       black_flag(CONCAT_MACRO(LEVEL,_black_flag_type),__FILE__,__LINE__,\
		  STRINGIFY(REVISION))) : 0)




#define WARN_DEBUG(REVISION,MESSAGE,LEVEL) \
(sprintf MESSAGE,\
 black_flag(CONCAT_MACRO(LEVEL,_black_flag_type),__FILE__,__LINE__,\
	    STRINGIFY(REVISION)))


#define WARN_MSG(LEVEL,MESSAGE,REVISION) \
(sprintf MESSAGE,\
 black_flag(CONCAT_MACRO(LEVEL,_black_flag_type),__FILE__,__LINE__,\
	    REVISION))



/*********************************************************************
  *
  *     OUT_OF_BOUNDS:
  *     checks if
  *       MIN <= INDEX < MAX
  *
  *     e.g.
  *     OUT_OF_BOUNDS(ivar,0,nvar)
  *
  ********************************************************************/

#define OUT_OF_BOUNDS(INDEX,MINIMUM_BOUND,MAXIMUM_BOUND) \
((INDEX)<(MINIMUM_BOUND) || (INDEX)>=(MAXIMUM_BOUND))

void handle_crash(int sigcode);
int handle_crash_init(void (*crash_fun)());
int black_flag(int bug_level,
	       char sourcefile[],
	       int sourceline,
	       char sourcefile_revision[]);

char *Program_name;
char *Program_version;

void black_flag_init(char progname[],char progversion[]);
void black_flag_init_args(int narg,char *arg[],char progversion[]);

#endif

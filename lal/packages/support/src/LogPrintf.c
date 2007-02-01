/*
 * Copyright (C) 2005 Reinhard Prix
 * 
 *  [partially based on the MSG_LOG class in BOINC:
 *  Copyright (C) 2005 University of California]
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the 
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
 *  MA  02111-1307  USA
 */

/**
 * \author Reinhard Prix
 * \date 2005
 * \file
 * \brief General-purpose log-message handling, controlled by logLevel independent of lalDebugLevel
 * 	  mostly modelled after the MSG_LOG class in BOINC.
 */

/* Windows version fixed by Bernd Machenschalk */

/*---------- INCLUDES ----------*/
#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include <time.h>

#include "LogPrintf.h"

/*---------- DEFINES ----------*/
#  ifndef __GNUC__
static volatile const char *name  = "$Id$";
#  else
static volatile const char __attribute__ ((unused)) *name  = "$Id$";
#  endif /* !__GNUC__ */


/* hardcoded for now, should be made configurable */
#define LogOutput 	stderr		

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- internal Global variables ----------*/



static int LogLevel = 0;	/* to be set with LogPrint_set_level() */

/*---------- internal prototypes ----------*/
static void LogPrintf_va (LogLevel_t level, const char* format, va_list va );

static const char * LogGetTimestamp (void);
static const char * LogTimeToString(double t);
static double LogDtime(void);

static const char *LogFormatLevel( LogLevel_t level );

/*==================== FUNCTION DEFINITIONS ====================*/

/* set the log-level to be used in this module.
 * (allow independence of lalDebugLevel!)
 */
void 
LogSetLevel(LogLevel_t level)
{
  LogLevel = level;
  return;
}

/** Verbatim output of the given log-message if LogPrintf_level >= level */
void
LogPrintfVerbatim (LogLevel_t level, const char* format, ...)
{
    va_list va;
    va_start(va, format);

    if ( LogLevel < level ) 
      return;

    /* simply print this to output  */
    vfprintf (LogOutput, format, va );
    fflush(LogOutput);

    va_end(va);

} /* LogPrintfVerbatim() */


/** prefix the log-message by a timestamp and level
 */
void
LogPrintf (LogLevel_t level, const char* format, ...)
{
  va_list va;
  va_start(va, format);

  LogPrintf_va ( level, format, va );

  va_end(va);

} /* LogPrintf() */


/** Low-level log-printing function: prefix message by timestamp if given. 
 */
void
LogPrintf_va (LogLevel_t level, const char* format, va_list va )
{
  if ( LogLevel < level ) 
    return;
     
  fprintf(LogOutput, "%s [%s]: ", LogGetTimestamp(), LogFormatLevel(level) );
  vfprintf(LogOutput, format, va);
  fflush(LogOutput);

  return;

} /* LogPrintf_va() */



/* taken from BOINC: return time-string for given unix-time
 */
const char * 
LogTimeToString ( double t ) 
{
  static char buf[100];
  char finer[16];
  time_t x = (time_t)t;
  struct tm* tm = localtime(&x);

  int hundreds_of_microseconds=(int)(10000*(t-(int)t));

  if (hundreds_of_microseconds == 10000) {
    /* paranoia -- this should never happen! */
    hundreds_of_microseconds=0;
    t+=1.0;
  }
  
  strftime(buf, sizeof(buf)-1, "%Y-%m-%d %H:%M:%S", tm);
  sprintf(finer, ".%04d", hundreds_of_microseconds);
  strcat(buf, finer);

  return buf;

} /* LogTimeToString() */


/** taken from BOINC's dtime():
 *  return time of day (seconds since 1970) as a double
 *
 * FIXME: windows-version certainly broken right now!
 */

#define EPOCHFILETIME_SEC (11644473600.)
#define TEN_MILLION 10000000.

static double 
LogDtime (void) 
{
#ifdef _MSC_VER

  /* Windows version of dtime() from BOINC.
     Compile switch is MS compiler macro,
     because I suspect gettimeofday should be present in MinGW */

  LARGE_INTEGER time;
  FILETIME sysTime;
  double t;
  GetSystemTimeAsFileTime(&sysTime);
  time.LowPart = sysTime.dwLowDateTime;
  time.HighPart = sysTime.dwHighDateTime;  // Time is in 100 ns units
  t = (double)time.QuadPart;    // Convert to 1 s units
  t /= TEN_MILLION;                /* In seconds */
  t -= EPOCHFILETIME_SEC;     /* Offset to the Epoch time */
  return t;

#else

  struct timeval tv;
  gettimeofday(&tv, 0);

  return tv.tv_sec + (tv.tv_usec/1.e6);

#endif
} /* LogDtime() */

/* returns static timestamps-string for 'now' */
static const char *
LogGetTimestamp (void)
{
  return ( LogTimeToString ( LogDtime() ) );

} /* LogGetTimestamp() */

/* returns static string characterizing this log-level */
static const char *
LogFormatLevel( LogLevel_t level )
{
  static char buf[100];

  switch ( level )
    {
    case LOG_CRITICAL:
      sprintf (buf, "CRITICAL" );
      break;
    case LOG_NORMAL:
      sprintf (buf, "normal");
      break;
    case LOG_DEBUG:
      sprintf (buf, "debug");
      break;
    case LOG_DETAIL:
      sprintf (buf, "detail");
      break;
    default:
      sprintf (buf, "unknown(%d)", level );
      break;
    } /* switch(level) */

  return buf;

} /* LogFormatLevel() */



/** Output gsl_matrix in octave-format, using the given format for the matrix-entries
 * return -1 on error, 0 if OK.
 */
int
XLALfprintfGSLmatrix ( FILE *fp, const char *fmt, const gsl_matrix *gij )
{
  int cols, rows;
  int i, j;

  /* check user input */
  if ( !gij || !fp || !fmt )
    return -1;

  rows = gij->size1;
  cols = gij->size2;

  fprintf (fp, " [ \\\n" );
  for ( i=0; i < rows; i ++ )
    {
      for (j=0; j < cols; j ++ )
	{
	  fprintf (fp, fmt, gsl_matrix_get ( gij, i, j ) );
	  if ( j < cols - 1 )
	    fprintf (fp, ", ");
	  else
	    fprintf (fp, ";\n");
	} /* for j < cols */
    }
  fprintf (fp, " ];\n" );
  
  return 0;

} /* XLALprintGSLmatrix() */


/** Output gsl_matrix in octave-format, using the given format for the matrix-entries
 * return -1 on error, 0 if OK.
 */
int
XLALfprintfGSLvector ( FILE *fp, const char *fmt, const gsl_vector *vect )
{
  int rows;
  int i;

  /* check user input */
  if ( !vect || !fp || !fmt )
    return -1;

  rows = vect->size;

  fprintf (fp, " [ " );
  for ( i=0; i < rows; i ++ )
    {
      fprintf (fp, fmt, gsl_vector_get ( vect, i ) );
      if ( i < rows - 1 )
	fprintf (fp, ", ");
    } /* for i < rows */

  fprintf (fp, " ];\n" );
  
  return 0;

} /* XLALprintGSLvector() */

int
XLALfprintfGSLvector_int ( FILE *fp, const char *fmt, const gsl_vector_int *vect )
{
  int rows;
  int i;

  /* check user input */
  if ( !vect || !fp || !fmt )
    return -1;

  rows = vect->size;

  fprintf (fp, " [ " );
  for ( i=0; i < rows; i ++ )
    {
      fprintf (fp, fmt, gsl_vector_int_get ( vect, i ) );
      if ( i < rows - 1 )
	fprintf (fp, ", ");
    } /* for i < rows */

  fprintf (fp, " ];\n" );
  
  return 0;

} /* XLALprintGSLvector_int() */

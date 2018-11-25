
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

int main( int argc, char **argv )
{
	FILE *ifile = stdin;
	char buf[256];
	if ( !fgets( buf, 256, ifile ) ) return 1;
	while ( fgets( buf, 256, ifile ) )
	{
		char *pp = strchr( buf, '.' );
		if ( !pp ) continue;
		*pp = 0;
		long year = atoi( buf );
		char *pp1 = strchr( pp + 1, '.' );
		*pp1 = 0;
		long month = atoi( pp + 1 );
		char *pp1d = strchr( pp1 + 1, ' ' );
		*pp1d = 0;
		long day = atoi( pp1 + 1 );
		char *pp1h = strchr( pp1d + 1, ':' );
		long hour = atoi( pp1d + 1 );
		char *pp1m = strchr( pp1h + 1, ':' );
		*pp1m = 0;
		long minute = atoi( pp1h + 1 );
		char *pp2 = strchr( pp1m + 1, ';' );
		double y = atof( pp2 + 1 );
		struct tm tms;
		tms.tm_sec = tms.tm_isdst = 0;
		tms.tm_min = minute;
		tms.tm_hour = hour;
		tms.tm_mday = day;
		tms.tm_mon = month;
		tms.tm_year = year - 1900;
		long x = mktime( &tms );
		printf( "%ld\t%g\n", x, y );
	}
	return 0;
}

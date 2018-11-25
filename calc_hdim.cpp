
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <string>
#include <set>
#include <map>
#include <vector>

using namespace std;

typedef pair<int,int> IPair;

const double big = 1e20;

struct DimEstimate
{
	double slope;
	double intercept;
	double correlation;
	int npoints;
	int dsize;
	int total;
	double residuals;
	double interval;
	vector<float> distr;
	DimEstimate() { npoints = dsize = total = 0; slope = intercept = correlation = residuals = interval = 0.; }
};

DimEstimate CalcDimension( vector< vector<float> >& distances, map<double,double>& tsums, map<double,int>& tcounts )
{
	vector<double> intervals;
	vector<double> radii;
	vector<double> rsums;
	int dsize = distances.size();
	double rsum = 0;
	for ( double interval = 0.005; interval < 0.05; interval += 0.001 )
	{
		int sump = 0;
		double vsump = 0.;
		for ( int c1 = 0; c1 < dsize; c1++ )
		{
			int numc = 1;
			double crsum = 0;
			for ( int c2 = 0; c2 < dsize; c2++ )
			{
				if ( c1 == c2 ) continue;
				if ( distances[ c1 ][ c2 ] < interval ) numc++;
			}
			sump += numc;
			vsump += log( dsize - numc );
		}
		if ( sump > 0 )
		{
			intervals.push_back( interval );
			double rad = double( sump ) / ( dsize * dsize );
			radii.push_back( rad );
			tsums[ interval ] += rad;
			tcounts[ interval ]++;
			rsums.push_back( sqrt( ( exp( vsump / dsize ) / dsize ) ) );
		}
	}
	DimEstimate rv;
	int npoints = intervals.size();
	if ( npoints < 3 ) return rv;
	double sumx = 0;
	double sumx2 = 0;
	double sumxy = 0;
	double sumy = 0;
	double sumy2 = 0;
	for ( int i = 0; i < npoints; i++ )
	{ 
		double x = intervals[i];
		double y = radii[i];
		sumx  += x;       
		sumx2 += x * x;  
		sumxy += x * y;
		sumy  += y;      
		sumy2 += y * y; 
		rv.distr.push_back( tsums[ intervals[i] ] / tcounts[ intervals[i] ] );
	} 
	double denom = ( npoints * sumx2 - sumx * sumx );
	rv.slope = ( npoints * sumxy  -  sumx * sumy ) / denom;
	rv.intercept = ( sumy * sumx2  -  sumx * sumxy ) / denom;
	rv.correlation = ( sumxy - sumx * sumy / npoints ) / sqrt( ( sumx2 - sumx * sumx/ npoints ) * ( sumy2 - sumy * sumy / npoints ) );
	rv.npoints = npoints;
	rv.dsize = dsize;
	double rlim = 0.9;
	double rinterval = ( rlim - rv.intercept ) / rv.slope;
	int iinterval = 1;
	while ( iinterval < intervals.size() && intervals[ iinterval ] < rinterval ) iinterval++;
	if ( iinterval >= 2 )
	{
		rv.residuals = rsums[ iinterval - 2 ] - rsums[ iinterval - 1 ];
		rv.interval = intervals[ iinterval - 1 ];
	}
	return rv;
}

bool CalcDistances( vector<string>& align, vector< vector<float> >& distances )
{
	distances.resize( align.size() );
	for ( int sc = 0; sc < align.size(); sc++ )
	{
		distances[sc].resize( align.size(), 1. );
		for ( int sc1 = 0; sc1 <= sc; sc1++ )
		{
			if ( align[sc1].size() != align[sc].size() )
			{
				fprintf( stderr, "sequences are not aligned\n" );
				return false;
			}
			int smuts = 0;
			int pg[2] = {0};
			int ml[2] = {0};
			for ( int pos = 0; pos < align[ sc ].size(); pos++ )
			{
				if ( align[sc][pos] == '-' && align[sc1][pos] =='-' ) continue;
				if ( align[sc][pos] == '-' )
				{
					pg[1] = 0;
					if ( pg[0] ) continue;
					else 
					{
						smuts++;
						pg[0] = 1;
						ml[1]++;
					}
				}
				else if ( align[sc1][pos] == '-' )
				{
					pg[0] = 0;
					if ( pg[1] ) continue;
					else 
					{
						smuts++;
						pg[1] = 1;
						ml[0]++;
					}
				}
				else
				{
					ml[0]++;
					ml[1]++;
					pg[0] = 0;
					pg[1] = 0;
					if ( toupper( align[sc][pos] ) != toupper( align[ sc1 ][ pos ] ) ) smuts++;
				}
			}
			int size = min( ml[0], ml[1] );
			double dist = double( smuts ) / size;
			distances[ sc ][ sc1 ] = dist;
			//fprintf( stderr, "%g ", dist );
		}
		//fprintf( stderr, "\n" );
	}
	for ( int sc = 0; sc < align.size(); sc++ )
	{
		for ( int sc1 = 0; sc1 < sc; sc1++ )
		{
			distances[ sc1 ][ sc ] = distances[ sc ][ sc1 ];
		}
	}
	return true;
}

const size_t alphabets = 26;

int min3( int a, int b, int c )
{
    if (a <= b && a <= c)
        return a;
    else if (b <= a && b <= c)
        return b;
    else
        return c;
}

int NWAlign(const string &a, const string &b, int alpha_gap, int alpha[alphabets][alphabets], string &a_aligned, string &b_aligned )
{
    size_t n = a.size();
    size_t m = b.size();

    vector<vector<int> > A(n + 1, vector<int>(m + 1));

    for (size_t i = 0; i <= m; ++i)
        A[0][i] = alpha_gap * i;
    for (size_t i = 0; i <= n; ++i)
        A[i][0] = alpha_gap * i;

    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];
            A[i][j] = min3( A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'],
                            A[i-1][j] + alpha_gap,
                            A[i][j-1] + alpha_gap );
        }
    }

    // print2DVector(A);

    a_aligned = "";
    b_aligned = "";
    size_t j = m;
    size_t i = n;
    for (; i >= 1 && j >= 1; )
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j] == A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'])
        {
            /*
             * I think prepending chars this way to a std::string is very inefficient.
             * Is there any better way of doing this without using C-style strings?
             */
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
			--i;
            --j;
        }
        else if (A[i][j] == A[i-1][j] + alpha_gap)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = '-' + b_aligned;
			--i;
        }
        else
        {
            a_aligned = '-' + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned = a[i-1] + a_aligned;
        b_aligned = '-' + b_aligned;
        --i;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned = '-' + a_aligned;
        b_aligned = b[j-1] + b_aligned;
        --j;
    }

    return A[n][m];
}

int CigarAlign(const string &rseq, const string &qseq, const string& cigar, int rbeg, string &r_aligned, string &q_aligned )
{
	r_aligned = rseq.substr( 0, rbeg );
	q_aligned = string( rbeg, '-' );
	char digits[16];
	char symb;
	int npos = 0;
	int ndigit = 0;
	int slen = cigar.size();
	int qbeg = 0;
	while ( npos < slen )
	{
		int symb = cigar[ npos ];
		if ( symb >= '0' && symb <= '9' )
		{
			digits[ ndigit++ ] = symb;
		}
		else
		{
			digits[ ndigit ] = 0;
			ndigit = 0;
			int val = atoi( digits );
			if ( symb == 'D' || symb == 'N' )
			{
				q_aligned += string( val, '-' );
				r_aligned += rseq.substr( rbeg, val );
				rbeg += val;
			}
			else if ( symb == 'I' )
			{
				q_aligned += qseq.substr( qbeg, val );
				r_aligned += string( val, '-' );
				qbeg += val;
			}
			else if ( symb == 'M' || symb == '=' || symb == 'X' )
			{
				q_aligned += qseq.substr( qbeg, val );
				r_aligned += rseq.substr( rbeg, val );
				qbeg += val;
				rbeg += val;
			}
		}
		npos++;
	}
	q_aligned += string( int( rseq.size() ) - rbeg, '-' );
	r_aligned += rseq.substr( rbeg );
}
	

bool LoadSeq( FILE *sfile, long pos, string& seq )
{
	fseek( sfile, pos, SEEK_SET );
	int state = 0;
	seq.clear();
	do
	{
		int symb = fgetc( sfile );
		if ( state == 0 )
		{
			if ( symb != '>' )
			{
				fprintf( stderr, "assert: '>' expected\n" );
				return false;
			}
			state = 1;
		}
		else if ( state == 1 )
		{
			if ( symb == EOF )
			{
				fprintf( stderr, "sequence expected at the end of file (pos=%ld ftell=%ld)\n", pos, ftell( sfile ) );
				return false;
			}
			if ( symb == '\n' )
			{
				state = 2;
			}
		}
		else if ( state == 2 )
		{
			if ( symb == EOF || symb == '>' ) break;
			if ( symb == '.' ) symb = '-';
			if ( isalpha( symb ) || symb == '-' ) seq.append( 1, toupper( symb ) );
		}
	}
	while ( 1 );
	return true;
}

bool GetHier( const string& id, vector<string>& nkey, bool eukflag )
{
	if ( id.find( ";" ) == -1 )
	{
		nkey.push_back( id );
		return true;
	}
	const int nlevel = 8;
	char tbuf[4096];
	strcpy( tbuf, id.data() );
	vector<string> key;
	key.resize( nlevel );
	const char *kw = "Bacteria;";
	const char *kw2 = "Archaea;";
	const char *kw3 = "Eukaryota;";
	const char *pkw = kw;
	char *pk = strstr( tbuf, kw );
	if ( !pk ) 
	{
		char *pk2 = strstr( tbuf, kw2 );
		if ( !pk2 ) 
		{
			char *pk3 = strstr( tbuf, kw3 );
			if ( !pk3 ) return false;
			else
			{
				pk = pk3;
				pkw = kw3;
				key[0] = kw3;
			}
		}
		else
		{
			pk = pk2;
			pkw = kw2;
			key[0] = "Archaea";
		}
	}
	else key[0] = "Bacteria";
	const char *pk0 = pk;
	pk += strlen( pkw );
	pk += strspn( pk, " \";" );
	char *cpk = pk;
	int lc = 1;
	int len0 = strlen( pk );
	for ( ; lc < nlevel; lc++ )
	{
		int nl = strcspn( cpk, "\";" );
		cpk[ nl ] = 0;
		int clc = lc;
		if ( strncmp( cpk, "otu", 3 ) == 0 ) break;
		key[ clc ] = cpk;
		if ( cpk + strlen( cpk ) + 1 - pk >= len0 ) break;
		char *npl = cpk + strlen( cpk ) + 1;
		npl += strspn( npl, " \";" );
		cpk = npl;
	}
	//nkey.resize( 5 );
	if ( eukflag ) 
	{
		nkey = key;
		int np = strcspn( id.data(), " " ) + 1;
		nkey.back() = string( id, np, pk0 - tbuf - np );
		return true;
	}
	
	if ( key[6].size() > 0 )
	{
		int kc;
		for (  kc = 3; kc < 6; kc++ ) if ( key[kc].substr( key[kc].size() - 3, 3 ) == "les" ) break;
		if ( kc == 6 ) nkey.insert( nkey.begin(), key.begin(), key.end() );
		else
		{
			nkey.insert( nkey.begin(), key.begin(), key.begin() + 2 );
			for ( int kcc = 3; kcc < kc; kcc++ ) nkey.back().append( ";" + key[kcc] );
			nkey.insert( nkey.end(), key.begin() + kc, key.end() );
		}
		while ( nkey.size() > 6 && nkey[6].size() )
		{
			nkey[4].append( ";" + nkey[5] );
			nkey.erase( nkey.begin() + 5 );
		}
		if ( nkey.size() > 6 ) nkey.erase( nkey.begin() + 6, nkey.end() );
	}
	else
	{
		nkey.insert( nkey.begin(), key.begin(), key.begin() + 6 );
	}
	return true;
}

void MergePairs( vector< vector<short> >& ralign, multimap<IPair,IPair>& qpairs, vector<string>& qseqs, vector<string>& qalignseqs )
{
	vector< vector<short> > qalign;
	qalign.resize( ralign.size() );
	int nqseqs = qseqs.size();
	int nrseqs = ralign[0].size();
	for ( int hc = 0; hc < ralign.size(); hc++ )
	{
		qalign[hc].resize( nqseqs, -1 );
		set<IPair> begs;
		for ( int bc = 0; bc < ralign[hc].size(); bc++ )
		{
			if ( ralign[hc][bc] != -1 ) begs.insert( IPair( bc, ralign[hc][bc] ) );
		}
		for ( set<IPair>::iterator it = begs.begin(); it != begs.end(); it++ )
		{
			for ( multimap<IPair,IPair>::iterator it1 = qpairs.lower_bound( *it ); it1 != qpairs.upper_bound( *it ); it1++ )
			{
				int qseqnum = it1->second.first - nrseqs;
				qalign[hc][ qseqnum ] = it1->second.second;
			}
		}
	}
	vector<int> begs( nqseqs, -1 );
	qalignseqs.resize( nqseqs );
	int bound = 0;
	for ( int hc = 0; hc < qalign.size(); hc++ )
	{
		int gmax = 0;
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			if ( qalign[hc][sc] != -1 )
			{
				gmax = max( gmax, bound + qalign[hc][sc] - begs[sc] );
			}
		}
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			if ( qalign[hc][sc] != -1 )
			{
				string sfragm( qseqs[sc], begs[sc] + 1, qalign[hc][sc] - begs[sc] );
				for ( int pc = 0; pc < int( sfragm.size() ) - 1; pc++ ) sfragm[pc] = tolower( sfragm[pc] );
				int gsize = ( gmax - qalignseqs[sc].size()  - sfragm.size() );
				sfragm.insert( 0, gsize, '-' );
				qalignseqs[sc] += sfragm;
				begs[sc] = qalign[hc][sc];
				bound = qalignseqs[sc].size();
			}
		}
		//for ( int sc = 0; sc < nqseqs; sc++ ) fprintf( stderr, "%3d %d %c ", begs[sc], qalign[hc][sc], ( qalignseqs[sc].size() == 0 ) ? '*' : qalignseqs[sc][ qalignseqs[sc].size() - 1 ] );
		//fprintf( stderr, " %d\n", hc );
	}
	{
		int gmax = 0;
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
				gmax = max( gmax, bound + int( qseqs[sc].size() ) - begs[sc] - 1 );
		}
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			string sfragm( qseqs[sc], begs[sc] + 1 );
			for ( int pc = 1; pc < int( sfragm.size() ); pc++ ) sfragm[pc] = tolower( sfragm[pc] );
			int gsize = ( gmax - qalignseqs[sc].size()  - sfragm.size() );
			sfragm.append( gsize, '-' );
			qalignseqs[sc] += sfragm;
		}
	}
}

struct RefData
{
	string ref_filename;
	vector<string> ids;
	vector<long> positions;
	map<string,int> index;
	map<int,string> otu;
	map<string,set<int> > rev_otu;
	bool LoadRefs( const char *filename, const char *kwstr );
};

bool RefData::LoadRefs( const char *filename, const char *kwstr )
{
	FILE *ifile = fopen( filename, "rt" );
	if ( !ifile )
	{
		fprintf( stderr, "can't open %s\n", filename );
		return false;
	}
	ref_filename = filename;
	char buf[65535];
	long fpos = 0;
	while ( fgets( buf, 65535, ifile ) )
	{
		strtok( buf, "\r\n" );
		if ( buf[0] == '>' )
		{
			if ( kwstr == 0 || strstr( buf, kwstr ) != 0 )
			{
				index[ string( buf + 1, 0, strcspn( buf + 1, " " ) ) ] = ids.size();
				char *potu = strstr( buf, "otu_" );
				string sotu;
				if ( potu )
				{
					sotu = string( potu, 0, strcspn( potu, " ;,.\n" ) );
				}
				else
				{
					if ( isdigit( buf[1] ) )
					{
						sotu = "otu_" + string( buf + 1, 0, strspn( buf + 1, "0123456789" ) );
					}
				}
				if ( sotu.size() )
				{
					otu[ ids.size() ] = sotu;
					rev_otu[ sotu ].insert( ids.size() );
				}
				ids.push_back( buf + 1 );
				positions.push_back( fpos );
			}
		}
		fpos = ftell( ifile );
	}
	fclose( ifile );
	fprintf( stderr, "input file loaded, index %d, otu %d, rotu %d\n", int( index.size() ), int( otu.size() ), int( rev_otu.size() ) );
	return true;
}

struct SamAlign
{
	string ref;
	int beg;
	string cigar;
};

struct SamData : RefData
{
	map<string,SamAlign> smap;
	map<string, map<int, set<string> > > otusets;
	set<string> rcreads;
	map<string,IPair> otubounds;
	bool LoadSam( const char *sname );
	bool LoadSamGen( const char *lname );
};

IPair parse_cigar( int ref_pos, const char *cigar, string& cleancigar )
{
	char digits[16];
	char symb;
	int npos = 0;
	int ndigit = 0;
	int slen = strlen( cigar );
	bool first = true;
	int beg = 0;
	int length = 0;
	int cbpos = 0;
	int cepos = 0;
	while ( npos < slen )
	{
		int symb = cigar[ npos ];
		if ( symb >= '0' && symb <= '9' )
		{
			digits[ ndigit++ ] = symb;
		}
		else
		{
			digits[ ndigit ] = 0;
			ndigit = 0;
			int val = atoi( digits );
			if ( symb == 'D' || symb == 'N' )
			{
				if ( first ) 
				{
					beg = val;
					cbpos = npos + 1;
				}
				else if ( npos + 1 == slen ) break;
				else 
				{
					length += val;
					cepos = npos;
				}
			}
			else if ( symb == 'M' || symb == '=' || symb == 'X' )
			{
				length += val;
				cepos = npos;
			}
			first = false;
		}
		npos++;
	}
	cleancigar = string( cigar + cbpos, 0, cepos - cbpos + 1 );
	return IPair( beg + ref_pos - 1, beg + length + ref_pos - 1 );
	
}

bool SamData::LoadSam( const char *sname )
{
	FILE *sfile = fopen( sname, "rt" );
	if ( !sfile ) 
	{
		fprintf( stderr, "can't open %s\n", sname );
		return false;
	}
	char buf[65535];
	while ( fgets( buf, 65535, sfile ) )
	{
		if ( buf[0] == '@' ) continue;
		vector<string> fields;
		for ( char *p = strtok( buf, "\t" ); p; p = strtok( 0, "\t" ) ) fields.push_back( p );
		string read_id( fields[0], 0, strcspn( fields[0].data(), " " ) );
		string ref_id( fields[2], 0, strcspn( fields[2].data(), " " ) );
		int flag = atoi( fields[1].data() );
		int ref_pos = atoi( fields[3].data() );
		string& cigar = fields[5];
		string cleancigar;
		IPair cbounds = parse_cigar( ref_pos, cigar.data(), cleancigar );
		if ( index.find( ref_id ) != index.end() ) 
		{
			int ref_ind = index[ ref_id ];
			string otuname = otu[ ref_ind ];
			SamAlign ca;
			ca.ref = ref_id;
			ca.beg = cbounds.first;
			ca.cigar = cleancigar;
			smap[ read_id ] = ca;
			otusets[ otuname ][ ref_ind ].insert( read_id );
			if ( flag & 16 ) rcreads.insert( read_id );
			if ( otubounds.find( otuname ) == otubounds.end() )
			{
				otubounds[ otuname ] = cbounds;
			}
			else
			{
				otubounds[ otuname ] = IPair( min( otubounds[ otuname ].first, cbounds.first ), max( otubounds[ otuname ].second, cbounds.second ) );
			}
		}
	}
	fclose( sfile );
	return true;
}

bool SamData::LoadSamGen( const char *lname )
{
	FILE *lfile = fopen( lname, "rt" );
	if ( !lfile ) 
	{
		fprintf( stderr, "can't open %s\n", lname );
		return false;
	}
	char buf[4096];
	bool first = true;
	while ( fgets( buf, 4096, lfile ) )
	{
		strtok( buf, "\r\n" );
		if ( strlen( buf ) == 0 || buf[0] == '\n' || buf[0] == '\r' ) continue;
		if ( first )
		{
			FILE *tfile = fopen( buf, "rt" );
			if ( !tfile )
			{
				LoadSam( lname );
				break;
			}
			fclose( tfile );
			first = false;
		}
		LoadSam( buf );
	}
	fclose( lfile );
	fprintf( stderr, "sam import done, %d motifs, %d otusets\n", int( smap.size() ), int( otusets.size() ) );
	fflush( stderr );
	return true;
}

struct ReadsData : SamData
{
	map<string,long> read_positions;
	map<string,string> read_files;
	map<string,int> samples_counts;
	bool LoadReads( const char *fname );
	bool LoadReadsGen( const char *fname );
	set<string> volumes;
};

string samplename( const string& fname )
{
	int pslash = fname.rfind( "/" );
	int pdot = fname.rfind( "." );
	if ( pslash != -1 )
	{
		if ( pdot != -1 && pdot > pslash ) return fname.substr( pslash + 1, pdot - pslash - 1 );
		return fname.substr( pslash + 1 );
	}
	else
	{
		if ( pdot != -1 ) return fname.substr( 0, pdot );
		return fname;
	}
}

bool ReadsData::LoadReads( const char *fname )
{
	FILE *ifile = fopen( fname, "rt" );
	if ( !ifile ) 
	{
		fprintf( stderr, "can't open %s\n", fname );
		return false;
	}
	volumes.insert( samplename( fname ) );
	char buf[65535];
	long fpos = 0;
	while ( fgets( buf, 65535, ifile ) )
	{
		strtok( buf, "\r\n" );
		if ( buf[0] == '>' )
		{
			string id( buf + 1, 0, strcspn( buf + 1, " \t" ) );
			if ( smap.find( id ) != smap.end() )
			{
				read_positions[ id ] = fpos;
				read_files[ id ] = fname;
				samples_counts[ samplename( fname ) ] ++;
			}
		}
		fpos = ftell( ifile );
	}
	fclose( ifile );
	return true;
}

bool ReadsData::LoadReadsGen( const char *lname )
{
	FILE *lfile = fopen( lname, "rt" );
	if ( !lfile ) 
	{
		fprintf( stderr, "can't open %s\n", lname );
		return false;
	}
	char buf[4096];
	bool first = true;
	while ( fgets( buf, 4096, lfile ) )
	{
		strtok( buf, "\r\n" );
		if ( strlen( buf ) == 0 || buf[0] == '\n' || buf[0] == '\r' ) continue;
		if ( first )
		{
			FILE *tfile = fopen( buf, "rt" );
			if ( !tfile )
			{
				LoadReads( lname );
				break;
			}
			fclose( tfile );
			first = false;
		}
		LoadReads( buf );
		fprintf( stderr, "file %s, %d /%d reads\n", buf, int( read_positions.size() ), int( read_files.size() )  );
	}
	fclose( lfile );
	fprintf( stderr, "reads import done, %d /%d reads\n", int( read_positions.size() ), int( read_files.size() )  );
	fflush( stderr );
	return true;
}
	
struct MainClass : ReadsData
{
	bool eukflag;
	map<int,int> cref_index;
	vector<int> cref_center;
	vector< map<int,string> > cref_seqs;
	vector< vector< vector<short> > > cref_aligns;
	FILE *of_reads;
	int clcnt1;
	int clcnt2;
	map< vector<string>, map<string,DimEstimate> > emap;
	
	MainClass() { clcnt1 = clcnt2 = 1; eukflag = false; }
	bool SingleRefs( const string& cotu );
	bool ProcessReads( const string& cotu, double threshold = 0.03, int maxgsize = 1000, const char *ofd_name = 0 );
};


bool MainClass::SingleRefs( const string& cotu )
{
	cref_index.clear();
	cref_center.clear();
	cref_seqs.clear();
	cref_aligns.clear();
	FILE *rfile = fopen( ref_filename.data(), "rt" );
	if ( !rfile )
	{
		fprintf( stderr, "assertion: rfile == 0 at group refs\n" );
		return false;
	}
	vector<string> rseqs;
	vector< set<int> > cl_index;
	for ( map<int,set<string> >::iterator it = otusets[ cotu ].begin(); it != otusets[ cotu ].end(); it++ )
	{
		rseqs.resize( rseqs.size() + 1 );
		if ( !LoadSeq( rfile, positions[ it->first ], rseqs.back() ) ) return false;
		set<int> singleset;
		singleset.insert( it->first );
		cl_index.push_back( singleset );
	}
	fclose( rfile );
	cref_center.resize( cl_index.size() );
	vector<int> ccenter( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		cref_center[ cc ] = *( cl_index[cc].begin() );
	}
	int bounds = 20;
	int rbeg = otubounds[ cotu ].first;
	int rend = otubounds[ cotu ].second;
	cref_seqs.resize( cl_index.size() );
	cref_aligns.resize( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		int abeg = rbeg;
		int aend = rend;
		for ( int lc = abeg; lc <= aend; lc++ )
		{
			vector<short> spos( 1, lc - abeg );
			cref_aligns[cc].push_back( spos );
		}
		cref_seqs[cc][cc] = rseqs[cc].substr( abeg, aend - abeg + 1 );
	}
	return true;
}

string itostr( int n )
{
	char buf[16];
	sprintf( buf, "%d", n );
	return string( buf );
}

bool MainClass::ProcessReads( const string& cotu, double threshold, int maxgsize, const char *ofd_name )
{
    int alpha[alphabets][alphabets];
    for (size_t i = 0; i < alphabets; ++i)
    {
        for (size_t j = 0; j < alphabets; ++j)
        {
            if (i == j) alpha[i][j] = 0;
            else alpha[i][j] = 1;
        }
    }
    int gap_penalty = 2;

	FILE *rfile = 0;
	string rfname;
	string sfname;
	int seqc = 0;
	map<double,double> tsums;
	map<double,int> tcounts;
	for ( map< int,set<string> >::iterator it0 = otusets[ cotu ].begin(); it0 != otusets[ cotu ].end(); it0++, seqc++ )
	{
		fprintf( stderr, "  ref %d, %d reads\n", seqc, int( it0->second.size() ) );
		vector< map<string,vector<string> > > rseqs;
		vector< map<string,multimap<IPair,IPair> > > ralign;
		vector< set<string> > rsamples;
		vector< vector<double> > rscores;
		vector< map<string,vector<string> > > rids;
		vector< map<string,vector< vector<short> > > > cralign;
		rseqs.resize( cref_seqs.size() );
		ralign.resize( cref_seqs.size() );
		rscores.resize( cref_seqs.size() );
		rsamples.resize( cref_seqs.size() );
		rids.resize( cref_seqs.size() );
		cralign.resize( cref_seqs.size() );
		for ( set<string>::iterator it = it0->second.begin(); it != it0->second.end(); it++ )
		{
			if ( read_files.find( *it ) == read_files.end() ) continue;
			if ( rfile == 0 || rfname != read_files[ *it ] )
			{
				if ( rfile ) fclose( rfile );
				rfname = read_files[ *it ];
				rfile = fopen( rfname.data(), "rt" );
				if ( !rfile )
				{
					fprintf( stderr, "can't open %s [%s]\n", rfname.data(), it->data() );
					return false;
				}
				sfname = samplename( rfname );
			}
			string rseq;
			if ( !LoadSeq( rfile, read_positions[ *it ], rseq ) ) return false;
			if ( rcreads.find( *it ) != rcreads.end() )
			{
				string rcseq;
				for ( int lc = 0; lc < rseq.size(); lc++ )
				{
					int clet = rseq[ rseq.size() - lc - 1 ];
					int nlet;
					switch ( clet )
					{
						case 'A': nlet = 'T'; break;
						case 'T': nlet = 'A'; break;
						case 'G': nlet = 'C'; break;
						case 'C': nlet = 'G'; break;
						nlet = 'N';
					}
					rcseq.append( 1, nlet );
				}
				rseq = rcseq;
			}
			int cluster = cref_index[ it0->first ];
			string qalign;
			string salign;
			string& sseq = cref_seqs[cluster][seqc];
			//int al_res = NWAlign( rseq, sseq, gap_penalty, alpha, qalign, salign );
			CigarAlign( sseq, rseq, smap[ *it ].cigar, smap[ *it ].beg - otubounds[ cotu ].first, salign, qalign );
			//fprintf( stderr, "%s\n%s\n\n", qalign.data(), salign.data() );
			int nseq = rseqs[ cluster ][ sfname ].size();
			rseqs[ cluster ][ sfname ].push_back( rseq );
			rids[ cluster ][ sfname ].push_back( *it );
			int b1 = 0;
			int b2 = 0;
			int nident = 0;
			int cref_size = cref_aligns[ cluster ][0].size(); 
			int cseqc = seqc;//cralign[ cluster ][ sfname ].size();
			for ( int lc = 0; lc < qalign.size(); lc++ )
			{
				if ( qalign[lc] == '-' ) b2++;
				else if ( salign[lc] == '-' ) b1++;
				else
				{
					IPair p1( cseqc, b2 );
					IPair p2( nseq + cref_size, b1 );
					ralign[ cluster ][ sfname ].insert( pair<IPair,IPair>( p1, p2 ) );
					ralign[ cluster ][ sfname ].insert( pair<IPair,IPair>( p2, p1 ) );
					if ( qalign[ lc ] == salign[ lc ] ) nident++;
					b1++;
					b2++;
				}
			}
			//fprintf( stderr, "%d %d\n\n", b1, b2 );
			rscores[ cluster ].push_back( double( nident ) / rseq.size() ); 
			rsamples[ cluster ].insert( sfname );
			//cralign[ cluster ][ sfname ].push_back( cref_aligns[ cluster ][seqc] );
			//if ( rseqs.size() >= maxgsize ) break;
		}
		for ( int cc = 0; cc < cref_seqs.size(); cc++ )
		{
			int otucl = clcnt1;
			clcnt1++;
			vector<string> hkey;
			GetHier( ids[ cref_center[cc] ], hkey, eukflag );
			hkey.push_back( itostr( otucl ) );
			for ( set<string>::iterator it = rsamples[ cc ].begin(); it != rsamples[ cc ].end(); it++ )
			{
				DimEstimate dimension;
				if ( rseqs[ cc ][ *it ].size() < maxgsize )
				{
					fprintf( stderr, "   %d reads to process\n", int( rseqs[ cc ][ *it ].size() ) );
					vector<string> ralignseqs;
					MergePairs( cref_aligns[cc], ralign[ cc ][ *it ], rseqs[ cc ][ *it ], ralignseqs );
					vector< vector<float> > distances;
					fprintf( stderr, "%s %s\n", hkey[0].data(), it->data() );
					if ( !CalcDistances( ralignseqs, distances ) ) return false;
					dimension = CalcDimension( distances, tsums, tcounts );
					dimension.total = samples_counts[ *it ];
					for ( int c2c = 0; c2c < ralignseqs.size(); c2c++ )
					{
						fprintf( of_reads, ">%s %s otu_%s\n%s\n", rids[cc][ *it ][c2c].data(), it->data(), hkey[0].data(), ralignseqs[c2c].data() );
					}
					fflush( of_reads );
				}
				emap[ hkey ][ *it ] = dimension;
			}
		}
	}
	if ( rfile ) fclose( rfile );
	if ( ofd_name )
	{
		FILE *odfile = fopen( ofd_name, "wt" );
		for ( map<double,double>::iterator it = tsums.begin(); it != tsums.end(); it++ ) 
		{
			fprintf( odfile, "%g %g\n", it->first, it->second / tcounts[ it->first ] );
		}
		fclose( odfile );
	}
	return true;
}

int main( int argc, char **argv )
{
	if ( argc == 1 )
	{
		fprintf( stderr, "arguments: ref_fasta sam/list reads/list out_table.txt out_align.fa\n" );
		return 1;
	}
	const char *refname = argv[1];
	const char *samname = argv[2];
	const char *readsname = argv[3];
	const char *of_name = argv[4];
	const char *ofa_name = argv[5];
	const char *kwstr = 0;
	double thresh1 = 0.6;
	double thresh2 = 0.3;
	int maxgsize = 50000;
	
	MainClass data;
	if ( !data.LoadRefs( refname, kwstr ) ) return 1;
	if ( !data.LoadSamGen( samname ) ) return 1;
	if ( !data.LoadReadsGen( readsname ) ) return 1;
	data.of_reads = fopen( ofa_name, "wt" );
	for ( map<string, map<int,set<string> > >::iterator it = data.otusets.begin(); it != data.otusets.end(); it++ )
	{
		fprintf( stderr, "%s: %d refs... ", it->first.data(), int( it->second.size() ) ); 
		fflush( stderr );
		//if ( !data.GroupRefs( it->first, thresh1, clmethod1 ) ) return 1;
		if ( !data.SingleRefs( it->first ) ) return 1;
		fprintf( stderr, "%d ref clusters\n", int( data.cref_seqs.size() ) );
		if ( !data.ProcessReads( it->first, thresh2, maxgsize, 0 ) ) return 1;
	}
	int nlevels = data.emap.begin()->first.size();
	FILE *of_table = fopen( of_name, "wt" );
	/*
	for ( int lc = 1; lc <= nlevels; lc++ ) fprintf( of_table, "%d\t", lc );
	for ( set<string>::iterator it = data.volumes.begin(); it != data.volumes.end(); it++ )
	{
		fprintf( of_table, "%s\t", it->data() );
	}
	fprintf( of_table, "\n" );
	*/
	for ( map< vector<string>, map<string,DimEstimate> >::iterator it = data.emap.begin(); it != data.emap.end(); it++ )
	{
		for ( set<string>::iterator it1 = data.volumes.begin(); it1 != data.volumes.end(); it1++ )
		{
			DimEstimate dim;
			if ( it->second.find( *it1 ) != it->second.end() ) dim = it->second[ *it1 ];
			for ( int lc = 0; lc < nlevels; lc++ ) fprintf( of_table, "%s\t", it->first[lc].data() );
			fprintf( of_table, "%s\t%g/%g/%g/%d/%d/%d/%g/%g\t", it1->data(), dim.slope, dim.intercept, dim.correlation, dim.npoints, dim.dsize, dim.total, dim.residuals, dim.interval );
			for ( int dc = 0; dc < dim.distr.size(); dc++ ) fprintf( of_table, "%g ", dim.distr[dc] );
			fprintf( of_table, "\n" );
		}
	}
	fclose( of_table );
	return 0;
}
	
	

	
	
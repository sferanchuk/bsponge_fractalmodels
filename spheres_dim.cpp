
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

enum { clSingle = 1, clComplete, clUpgma };

typedef pair<int,int> IPair;

const double big = 1e20;


bool CalcDistances( map<string,string>& align, map< string,map<string,float> >& distances )
{
	for ( map<string,string>::iterator it1 = align.begin(); it1 != align.end(); it1++ )
	{
		for ( map<string,string>::iterator it2 = align.begin(); it2 != it1; it2++ )
		{
			if ( it1->second.size() != it2->second.size() )
			{
				fprintf( stderr, "sequences are not aligned\n" );
				return false;
			}
			int smuts = 0;
			int pg[2] = {0};
			int ml[2] = {0};
			for ( int pos = 0; pos < it1->second.size(); pos++ )
			{
				if ( it1->second[pos] == '-' && it2->second[pos] =='-' ) continue;
				if ( it1->second[pos] == '-' )
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
				else if ( it2->second[pos] == '-' )
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
					if ( toupper( it1->second[pos] ) != toupper( it2->second[ pos ] ) ) smuts++;
				}
			}
			int size = min( ml[0], ml[1] );
			double dist = double( smuts ) / size;
			distances[ it1->first ][ it2->first ] = dist;
			distances[ it2->first ][ it1->first ] = dist;
			fprintf( stderr, "[%d %d %d %d %g] ", smuts, ml[0], ml[1], int( it1->second.size() ), dist );
		}
		fprintf( stderr, "\n" );
	}
	return true;
}

void CalcTHist( map<string,set<string> >& iclusters, map<string, string>& seqs, map< string, map< string, float > >& dist, vector<int>& thist )
{
	vector< set<string> > cclusters;
	set<string> ids;
	for ( map<string,set<string> >::iterator it = iclusters.begin(); it != iclusters.end(); it++ ) 
	{
		ids.insert( it->second.begin(), it->second.end() );
		if ( it->second.size() > 2 ) cclusters.push_back( it->second );
	}
	const int hsize = 200;
	map<int,set<string> > clmap;
	for ( set<string>::iterator it = ids.begin(); it != ids.end(); it++ )
	{
		if ( dist.find( *it ) == dist.end() ) continue;
		double mindist = 100;
		int cclust = -1;
		for ( int cc = 0; cc < cclusters.size(); cc++ )
		{
			int dcnt = 0;
			double dsum = 0;
			for ( set<string>::iterator it1 = cclusters[cc].begin(); it1 != cclusters[cc].end(); it1++ )
			{
				if ( dist[ *it ].find( *it1 ) != dist[ *it ].end() )
				{
					dsum += dist[ *it ][ *it1 ];
					dcnt++;
				}
			}
			if ( dcnt > 0 )
			{
				double cdist = dsum / dcnt;
				if ( cdist < mindist ) 
				{
					mindist = cdist;
					cclust = cc;
				}
			}
		}
		if ( cclust != -1 ) clmap[ cclust ].insert( *it );
	}
	/*
	//map<string,string> seqs;
	//map<string,pair<int,int> > seqpos;
		if ( pseq.size() )
		{
			strcpy( buf, pseq.data() );
			int blen = strlen( buf );
			int beg = 0; 
			for ( ; ( buf[beg] == '-' || buf[beg] == '.' ) && beg < blen; beg++ );
			int end = blen - 1;
			for ( ; ( buf[end] == '-' || buf[end] == '.' ) && end >= 0; end-- );
			if ( end > beg )
			{
				seqs[cid] = buf;
				seqpos[cid] = pair<int,int>( beg, end );
			}
		}
		cid = nid;
		pseq.clear();
		}
		else
		{
			strtok( buf, "\n" );
			pseq.append( buf );
		}
	}
	if ( pseq.size() )
	{
		strcpy( buf, pseq.data() );
		int blen = strlen( buf );
		int beg = 0; 
		for ( ; ( buf[beg] == '-' || buf[beg] == '.' ) && beg < blen; beg++ );
		int end = blen - 1;
		for ( ; ( buf[end] == '-' || buf[end] == '.' ) && end >= 0; end-- );
		if ( end > beg )
		{
			seqs[cid] = buf;
			seqpos[cid] = pair<int,int>( beg, end );
		}
	}
	fclose( afile );
	*/
	vector<string> ideal;
	vector<int> cldim;
	ideal.resize( cclusters.size() );
	cldim.resize( cclusters.size(), 0 );
	thist.resize( hsize, 0 );
	for ( int cc = 0; cc < cclusters.size(); cc++ )
	{
		int slen = 0;
		map<int, map<int,int> > count;
		for ( set<string>::iterator it = cclusters[cc].begin(); it != cclusters[cc].end(); it++ )
		{
			if ( seqs.find( *it ) == seqs.end() ) continue;
			if ( slen == 0 ) slen = seqs[ *it ].size();
			//for ( int lc = seqpos[ *it ].first; lc <= seqpos[ *it ].second; lc++ ) count[ lc ][ seqs[ *it ][ lc ] ]++;
			for ( int lc = 0; lc <= seqs[ *it ].size(); lc++ ) count[ lc ][ seqs[ *it ][ lc ] ]++;
		}
		string cideal( slen, '-' );
		int cdim = 0;
		for ( map<int, map<int,int> >::iterator it = count.begin(); it != count.end(); it++ )
		{
			if ( it->second.size() == 1 ) 
			{
				cideal[ it->first ] = it->second.begin()->first;
			}
			else
			{
				cdim++;
				int cmax = 0;
				for ( map<int,int>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++ )
				{
					if ( it1->second > cmax )
					{
						cmax = it1->second;
						cideal[ it->first ] = it1->first;
					}
				}
			}
		}
		ideal[cc] = cideal;
		cldim[cc] = cdim;
		if ( clmap.find( cc ) == clmap.end() ) continue;
		int hist[hsize] = { 0 };
		for ( set<string>::iterator it = clmap[cc].begin(); it != clmap[cc].end(); it++ )
		{
			if ( seqs.find( *it ) == seqs.end() ) continue;
			int nmut = 0;
			for ( int lc = 0; lc <= seqs[ *it ].size(); lc++ ) 
			{
				if ( seqs[ *it ][ lc ] != cideal[lc] ) nmut++;
			}
			if ( nmut < hsize ) hist[ nmut ]++;
		}
		double lsumx = 0;
		double lsumxx = 0;
		double lsumy = 0;
		double lsumxy = 0;
		int lcnt = 0;
		for ( int i = 0; i < hsize; i++ )
		{
			//printf( "cl %d %d %d\n", cc, i, hist[i] );
			thist[i] += hist[i];
			if ( hist[i] > 0 )
			{
				double x = i;
				double y = log( hist[i] );
				lsumx += x;
				lsumy += y;
				lsumxx += x * x;
				lsumxy += x * y;
				lcnt++;
			}
		}
		if ( lcnt < 3 ) continue;
		double mx = lsumx / lcnt;
		double vdim = ( lsumxy - lsumy * mx )  / ( lsumxx - lsumx * mx );
		double y0 = lsumy / lcnt - vdim * lsumx / lcnt;
		double sumsq = 0;
		for ( int i = 0; i < hsize; i++ )
		{
			if ( hist[i] > 0 )
			{
				double x = i;
				double y = log( hist[i] );
				double ym = y0 + vdim * x;
				sumsq += ( ym - y ) * ( ym - y );
			}
		}
		double disp = sqrt( sumsq / lcnt ); 

//		printf( "sp %d %d %g %d %d %g\n", cc, cdim, vdim, int( cclusters[cc].size() ), int( clmap[cc].size() ), disp );
	}
}
	/*
	for ( int i = 0; i < hsize; i++ )
	{
		printf( "cl all %d %d\n", i, thist[i] );
	}
	{
		double lsumx = 0;
		double lsumxx = 0;
		double lsumy = 0;
		double lsumxy = 0;
		int lcnt = 0;
		for ( int i = 20; i < 100; i++ )
		{
			if ( thist[i] > 0 )
			{
				double x = i;
				double y = log( thist[i] );
				lsumx += x;
				lsumy += y;
				lsumxx += x * x;
				lsumxy += x * y;
				lcnt++;
			}
		}
		double mx = lsumx / lcnt;
		double vdim = ( lsumxy - lsumy * mx )  / ( lsumxx - lsumx * mx );
		double y0 = lsumy / lcnt - vdim * lsumx / lcnt;
		double sumsq = 0;
		for ( int i = 20; i < 100; i++ )
		{
			if ( thist[i] > 0 )
			{
				double x = i;
				double y = log( thist[i] );
				double ym = y0 + vdim * x;
				sumsq += ( ym - y ) * ( ym - y );
			}
		}
		double disp = sqrt( sumsq / lcnt ); 
		printf( "res dim %g %g\n", vdim, disp );
	}
	*/

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
	DimEstimate() { npoints = dsize = total = 0; slope = intercept = correlation = residuals = interval = 0.; }
};

DimEstimate CalcDimension( map<string,map<string,float> >& distances, map<double,double>& tsums, map<double,int>& tcounts )
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
		for ( map<string,map<string,float> >::iterator it1 = distances.begin(); it1 != distances.end(); it1++ )
		{
			int numc = 1;
			double crsum = 0;
			for ( map<string,float>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++ )
			{
				if ( it1->first == it2->first ) continue;
				if ( it2->second < interval ) numc++;
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

int main( int argc, char **argv )
{
	if ( argc < 2 )
	{
		fprintf( stderr, "arguments: alignment [-bysample]\n" );
		return 1;
	}
	FILE *ifile = fopen( argv[1], "rt" );
	int option = ( argc > 2 && strcmp( argv[2], "-bysample" ) == 0 ) ? 1 : 0;
	map<string, map<string, set<string> > > gids;
	map<string, map<string,string> > gseqs;
	set<string> refs;
	map<string,int> scounts;
	char buf[ 65535 ];
	string cid;
	string ccluster;
	string cref;
	while ( fgets( buf, 65535, ifile ) )
	{
		if ( buf[0] == '>' )
		{
			cid = string( buf + 1, 0, strcspn( buf + 1, " " ) );
			string sample = string( buf + 2 + cid.size(), 0, strcspn( buf + 2 + cid.size(), " " ) ); 
			char *potu = strstr( buf, "otu_" );
			string cotu = string( potu + 4, 0, strcspn( potu + 4, " \n" ) );
			refs.insert( cotu );
			cref = cotu;
			if ( option == 1 )
			{
				char nbuf[64];
				sprintf( nbuf, "%s\t%d\t%s", cotu.data(), int( refs.size() ), sample.data() );
				cref = nbuf;
			}
			//if ( gseqs[ cref ].size() > 2000 ) break;
			char *pcluster = strchr( potu + 4, '_' );
			if ( pcluster != 0 ) ccluster = string( pcluster + 1, 0, strcspn( pcluster + 1, "_" ) );
			gids[ cref ][ ccluster ].insert( cid );
			scounts[ sample ] ++;
		}
		else
		{
			strtok( buf, "\n" );
			gseqs[ cref ][ cid ].append( buf );
		}
	}
	for ( map<string, map<string, set<string> > >::iterator it0 = gids.begin(); it0 != gids.end(); it0++ )
	{
		map<string,set<string> >& iclusters = it0->second;
		map<string,string>& seqs = gseqs[ it0->first ];
		map< string, map< string, float > > dist;
		CalcDistances( seqs, dist );
		if ( option != 1 )
		{
			vector<int> thist;
			CalcTHist( iclusters, seqs, dist, thist );
			printf( "%s\tsdim\t", it0->first.data() );
			for ( int hc = 0; hc < thist.size(); hc++ ) printf( "%d ", thist[ hc ] );
			printf( "\n" );
		}
		map<double,double> tsums;
		map<double,int> tcounts;
		DimEstimate dim = CalcDimension( dist, tsums, tcounts );
		if ( option != 1 )
		{
			printf( "%s\t%g\t", it0->first.data(), dim.residuals );
		}
		else
		{
			string sample = it0->first.substr( it0->first.find( "\t" ) + 3 );
			printf( "%s\t%g/%g/%g/%d/%d/%d/%g/%g\t", it0->first.data(), dim.slope, dim.intercept, dim.correlation, dim.npoints, dim.dsize, scounts[ sample ], dim.residuals, dim.interval );
		}
		for ( map<double,double>::iterator it = tsums.begin(); it != tsums.end(); it++ ) printf( "%g ", it->second / tcounts[ it->first ] );
		printf( "\n" );
	}
	return 0;
}

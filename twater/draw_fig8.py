
import sys
from dateutil.parser import parse
from datetime import datetime
import time
import pyrem
import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
import math

sys.path.append( '/home/sergey/soft/HiguchiFractalDimension' )
import hfd


tinterval = 120
maxheat = 23.

gs = gridspec.GridSpec(5, 2, height_ratios=[ 4, 1, 4, 1, 0.5 ], width_ratios = [ 1, 1 ] )
fig = plt.figure( figsize=(20,20),dpi=80 ) 
#ax0 = plt.subplot( gs[0] )
#ax1 = plt.subplot( gs[1] )

#fig, ax = plt.subplots( 2, 4, figsize=(40,20), dpi=80 )	

locations = [ "listv", "bkot", "baik", "uzur" ]
locnames = [ "Listvyanka", "B.Koty", "Baikalsk", "Uzur" ]

xlpos = []
xlabels = []

for year in range( 2010, 2018 ):
	dt = datetime( year=year, month=1, day=1 )
	tstamp = time.mktime( dt.timetuple() )
	xlpos.append( tstamp )
	xlabels.append( str( year ) )
	
def moving_average( tseries, N ):
	res = np.convolve( tseries, np.ones((N,))/N, mode='valid' )
	return res


for lnum in range( len( locations ) ):
	idata = []
	itime = []
	with open( "twater_" + locations[ lnum ] + "_dig.txt" ) as f:
		for line in f:
			ls = line.split()
			itime.append( int( ls[0] ) )
			idata.append( float( ls[1] ) )
			continue
			dtstr = line[:19]
			tmstr = line[20:25]
			if dtstr[0] == 'D':
				continue
			#print dtstr
			dtime = parse( dtstr )
			#print dtime
			tstamp = time.mktime( dtime.timetuple() )
			itime.append( tstamp )
			idata.append( float( tmstr ) )
			#print ( tstamp, float( tmstr ) )
			#if len( idata ) > 200000:
			#	break


	hdata = []
	htime = []
	irdata = []
	hsize = 8000
	stepsize = 500
	idmax = -100.
	idmin = 100.
	minssize = 50

	chdata = []
	chtime = []
	cidata = []
	sumerr = 0.
	sumdim = 0
	sumcsize = 0
	minpv = [ 1 ] * 13
	maxpv = [ 0 ] * 13
	sumpv = [ 0. ] * 13
	sumpvc = [ 0 ] * 13

	for pnum in range( 0, len( idata ) - hsize, stepsize ):
		dsegment = idata[ pnum : pnum + hsize ]
		davg = np.average( dsegment )
		dstd = np.std( dsegment )
		tbreak = False
		for tk in range( pnum, pnum + hsize ):
			if tk > pnum and abs( itime[ tk - 1 ] - itime[ tk ] ) > 60 * 24 * 60 * 10:
				tbreak = True
			if abs( idata[ tk ] - davg ) > 10 * dstd:
				dsegment[ tk - pnum ] = davg
		if tbreak == True:
			if len( cidata ) > 0:
				hdata.append( chdata )
				htime.append( chtime )
				irdata.append( cidata )
				idmax = max( idmax, max( cidata ) )
				idmin = min( idmin, min( cidata ) )
			chdata = []
			chtime = []
			cidata = []
			continue
		chtime.append( itime[ pnum ] )
		( dim, err, csize, pv ) = hfd.hfd_ext( dsegment )
		cdt = datetime.utcfromtimestamp( itime[ pnum ] )
		minpv[ cdt.month ] = min( minpv[ cdt.month ], pv )
		maxpv[ cdt.month ] = max( maxpv[ cdt.month ], pv )
		if pv > 0 and pv < 1:
			sumpv[ cdt.month ] += math.log( pv )
			sumpvc[ cdt.month ] += 1
		
		chdata.append( dim )
		if not np.isnan( err ):
			sumerr += err
			sumdim += 1
			sumcsize += csize
		cidata.append( np.average( dsegment ) )

	hdata.append( chdata )
	htime.append( chtime )
	irdata.append( cidata )

	print >> sys.stderr, ( len( htime ), sumdim, sumerr, sumerr / sumdim, sumcsize / sumdim )
	print >> sys.stderr, minpv
	print >> sys.stderr, maxpv
	print >> sys.stderr, sumpv
	print >> sys.stderr, sumpvc
	xmin = min( xlpos )
	xmax = max( xlpos )
	
	ax1 = plt.subplot( gs[ ( 2 * ( lnum // 2 ) + 1, lnum % 2 ) ] )
	ax2 = plt.subplot( gs[ ( 2 * ( lnum // 2 ), lnum % 2 ) ] )
	ntidata = np.array( [ -1. ] * int( ( xmax - xmin ) / ( tinterval * stepsize ) ) )
	for snum in range( len( htime ) ):
		#ax1.plot( htime[ snum ], irdata[ snum ], color="steelblue" )
		for ct in range( len( htime[ snum ] ) ):
			ind = int( ( htime[ snum ][ ct ] - xmin ) / ( tinterval * stepsize ) )
			ntidata[ ind ] = irdata[ snum ][ ct ] / maxheat
			#print ind
			ax2.plot( htime[ snum ], hdata[ snum ], color="lightblue" )
			if len( htime[snum ] ) > minssize:
				ma = moving_average( hdata[ snum ], minssize )
				ax2.plot( htime[ snum ][ minssize / 2 : len( ma ) + minssize / 2 ], ma, color="steelblue" )
	print len( ntidata )
	print ntidata
	#mtidata = np.ma.masked_less( ntidata, 0 )
	#ccolmap = plt.get_cmap( "jet" )
	ccolmap = copy( plt.get_cmap( "jet" ) )
	#print mtidata
	ccolmap.set_bad( alpha = 0.0 )
	ccolmap.set_under( alpha = 0.0  )
	cs = ax1.imshow( np.vstack( ( ntidata, ntidata ) ), cmap = ccolmap, aspect = 0.015 * len( ntidata ), vmin = 0, vmax = 1 )# 'auto' )
	ax1.set_xticks( [] )
	#ax1.set_xticklabels( xlabels )
	ax2.set_xticks( xlpos )
	ax2.set_xticklabels( xlabels )
	for xp in xlpos:
		ax2.axvline( xp, color = "darkgrey" )
	ax2.axhline( 1.5, color = "black" )
	ax2.axhline( 1.3, color = "darkgreen" )
	#ax1.set_xlim( [ xmin, xmax ] )
	ax2.set_xlim( [ xmin, xmax ] )
	#ax1.set_ylim( [ 0, 22 ] )
	ax1.set_yticks( [] )
	ax2.set_ylim( [ 0.9, 2.1 ] )
	ax2.set_yticks( [ 1.1, 1.3, 1.5, 1.7, 1.9 ] )
	#ax1.set_title( locations[ lnum ] + " - temperature" )
	ax2.set_title( locnames[ lnum ], size = 24., loc="left" )# + " - fractal dimension (Higuchi method)" )
	ax2.tick_params( axis='y', labelsize = 24. )
	ax2.tick_params( axis='x', labelsize = 18. )

cbaxis = plt.subplot( gs[ (4,0) ] )
cbar = fig.colorbar( cs, cax = cbaxis, orientation = "horizontal", ticks = [ 0, 4./ maxheat, 10. / maxheat, 20. / maxheat ] )
cbar.ax.set_xticklabels( [ "0", "4", "10", "20" ] )
cbar.ax.tick_params( axis='x', labelsize = 24. )
cbar.ax.set_title( "Temperature", size = 24., loc="left" )
	
plt.savefig( "fig8.png" )

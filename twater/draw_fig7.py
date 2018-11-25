
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

sys.path.append( '/home/sergey/soft/HiguchiFractalDimension' )
import hfd


tinterval = 120
maxheat = 23.

gs = gridspec.GridSpec( 1, 2, width_ratios=[ 2 ] * 2 )
fig = plt.figure( figsize=(20,10),dpi=80 ) 
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

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


for lnum in range( 1 ):
	idata = []
	itime = []
	with open( "twater_" + locations[ lnum ] + "_dig.txt" ) as f:
		for line in f:
			ls = line.split()
			itime.append( int( ls[0] ) )
			idata.append( float( ls[1] ) )
			if len( idata ) > 100000:
				break
	for beg in [ 0, 50000 ]:
		ax = ax1 if beg == 0 else ax2
		bsize = 20000
		for wsize in [ 360, 720, 1440 ]:
			ma = moving_average( idata[ beg : beg + bsize ], wsize )
			psize = bsize / wsize
			ax.plot( [ float( r ) / psize for r in range( psize ) ], [ ma[ r * wsize ] for r in range( psize ) ], label = str( wsize / 30 ) + " h"  )
			xticks = [ 0, 0.5, 1 ]
		ax.set_xticks( xticks )
		ax.set_xticklabels( [ datetime.utcfromtimestamp( itime[ beg + int( xt * bsize ) ] ).strftime( '%d.%m.%Y' ) for xt in xticks ] )
		ax.tick_params( axis='y', labelsize = 20. )
		ax.tick_params( axis='x', labelsize = 20. )
ax1.set_ylabel( "Temperature", size = 20. )
ax1.legend( loc = "upper right", fontsize = 20 )
plt.savefig( "fig7.png" )

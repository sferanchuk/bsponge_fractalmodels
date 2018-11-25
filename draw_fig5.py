


import sys
from dateutil.parser import parse
from datetime import datetime
import time
import pyrem
import matplotlib
from matplotlib import gridspec
from matplotlib import colors as mcolors
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
import math
import scipy.stats as stats

snames = { "1118847" : "Trebouxiophyceae", "840888" : "Sediminibacterium", "848608" : "Rhodospirillaceae", "550168" : "Synechococcus", "4324048" : "Flavobacterium", "72607" : "Comamonadaceae", "217320" : "Rhodobacter" }
snkeys = snames.keys()

gs = gridspec.GridSpec( 1, 2, width_ratios=[ 2 ] * 2 )
fig = plt.figure( figsize=(20,15),dpi=80 ) 
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

sdistr = [ ( [], [] ) ] * len( snkeys )
spec = [ 0. ] * len( snkeys )
scolors = [ "steelblue", "indianred", "olive", "burlywood", "plum", "pink", "slategrey" ]

max_ident = 0.05
mid_ident = 0.02
min_ident = 0.003

def calc_jakovenko( sndata ):
	distr = sndata
	dmaxx = len( sndata )
	datax = []
	datay = []
	xn = float( len( distr ) )
	yn = float( distr[0] )
	for i in range( len( distr ) ):
		vy = math.log( distr[i] )
		vx = math.log( i + 1 )
		datax.append( vx )
		datay.append( vy )
	nmin = 3
	nmax = int( len( distr ) * 0.7 )
	rbest = 0
	ebest = 1000
	nbest = nmax
	pvbest = 1
	for cnsegm in range( nmin, nmax ):
		(ap_s,bp_s,rp,ttp,stderrp) = stats.linregress( datax[:cnsegm], datay[:cnsegm] )
		(ab_s,bb_s,rb,ttb,stderrb) = stats.linregress( range( cnsegm + 1, len( distr ) + 1 ), datay[ cnsegm : ] )
		yres = [ ap_s * math.log( mx ) + bp_s for mx in range( 1, cnsegm + 1 ) ] + [ ab_s * mx + bb_s for mx in range( cnsegm + 1, len( distr ) + 1 ) ]
		etot = np.std( [ yres[i] - datay[i] for i in range( len( distr ) ) ] )
		pvtot = 1 - ( 1 - ttp ) * ( 1 - ttb )
		rtot = rp + rb
		#etot = ( stderrp * math.sqrt( cnsegm ) + stderrb * math.sqrt( len( distr ) + 1 - cnsegm ) ) / math.sqrt( len( distr ) )
		#etot = ( stderrp * cnsegm + stderrb * ( len( distr ) + 1 - cnsegm ) ) / len( distr )
		if etot < ebest:
			pvbest = pvtot
			rbest = rtot
			ebest = etot
			ebestp = ( stderrp, stderrb, etot, cnsegm, rp, rb )
			nbest = cnsegm
			abestp = ap_s
			bbestp = bp_s
			abestb = ab_s
			bbestb = bb_s
			ybest = yres
	print ( abestp, bbestp, abestb, bbestb, pvbest, rbest, ebestp, len( sndata ), abestb * len( sndata ) / ( max_ident - min_ident ), 100. * ( min_ident + float( ebestp[3] ) * ( max_ident - min_ident ) / len( sndata ) ) )
	return ( range( 1, len( distr ) +  1 ), ybest, abestb, abestp )


yshift = 1.5

with open( "hdim_sph_table.txt" ) as f:
	for line in f:
		ls = line.strip().split( "\t" )
		snind = snkeys.index( ls[0] )
		snshift = snind * yshift
		sndata = map( float, ls[2].split() )
		if ls[1] != "sdim":
			print snames[ ls[0] ]
			( vx, vy, abest, abestl ) = calc_jakovenko( sndata )
			ax2.plot( vx, [ cy + snshift for cy in vy ], c = "black", lw = 0.7 ) #, c = scolors[ snind ] )
			ax2.scatter( range( 1, len( sndata ) + 1 ), [ math.log( v ) + snshift for v in sndata ], c = scolors[ snind ] )
			#ax1.scatter( range( 1, len( sndata ) + 1 ), [ math.log( v ) for v in sndata ], c = scolors[ snind ], label = "%s %5.3f %5.3f" % ( snames[ ls[0] ], abest, abestl ) )
			nvx = [ math.log( v ) for v in range( 1, len( sndata ) + 1 ) ]
			ax1.scatter( nvx, [ math.log( v ) + snshift for v in sndata ], c = scolors[ snind ], label = snames[ ls[0] ] )
			ax1.plot( nvx, [ cy + snshift for cy in vy ], c = "black", lw = 0.7 ) #scolors[ snind ] )
			ax2.text( len( sndata ), snshift + 0.1, snames[ ls[0] ], ha = "right", va = "bottom", size = 16 )
			
xlen = len( sndata )
#ax1.legend( loc = "lower right", fontsize = "16" )
ax1.set_ylim( [ -5, yshift * len( snkeys ) ] )
ax2.set_ylim( [ -5, yshift * len( snkeys ) ] )
ax1.set_yticks( [] )
ax2.set_yticks( [] )
#ax1.set_yticks( [ math.log( 0.01 ), math.log( 0.1 ), 0 ] )
#ax2.set_yticklabels( [ "1%", "10%", "100%" ] )
#ax2.set_yticks( [ math.log( 0.01 ), math.log( 0.1 ), 0 ] )
#ax1.set_yticklabels( [ "1%", "10%", "100%" ] )
ax1.set_xticks( [ 0, math.log( xlen * mid_ident / max_ident ), math.log( xlen ) ] )
ax1.set_xticklabels( [ str( min_ident ), str( mid_ident ), str( max_ident ) ] )
ax2.set_xticks( [ 0, xlen * mid_ident / max_ident, xlen ] )
ax2.set_xticklabels( [ str( min_ident ), str( mid_ident ), str( max_ident ) ] )
ax1.set_xlabel( "Distance (logarithmic scale)", size = 20. )
ax2.set_xlabel( "Distance (linear scale)", size = 20. )
ax1.set_ylabel( "Count of distances (logarithmic scale)", size = 20. )
ax1.tick_params( axis='y', labelsize = 16. )
ax1.tick_params( axis='x', labelsize = 16. )
ax2.tick_params( axis='y', labelsize = 16. )
ax2.tick_params( axis='x', labelsize = 16. )
plt.savefig( "fig5.png" )
			


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
snkeys = [ "1118847", "840888", "848608", "550168", "4324048", "72607", "217320" ]
pscales = [ 10.7, 4.9, 7.32, 8.2, 22.0, 28.3, 10.7 ]

shealthy = [ "L1", "L2", "L3", "L7", "L8", "OV1", "OV3", "OV4", "OV6", "T1", "T4", "T6" ]
llet = [ "L", "O", "T" ]

gs = gridspec.GridSpec( 2, 4, height_ratios=[ 2 ] * 2, width_ratios = [ 2 ] * 4 )
fig = plt.figure( figsize=(20,20),dpi=80 ) 

mmax = 0.05
#mmin = 1.

colmap = copy( plt.get_cmap( "copper" ) )
#colmap.set_under( "darkgreen", alpha = 0.0  )


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
	#print ( abestp, bbestp, abestb, bbestb, pvbest, rbest, ebestp )
	return ( range( 1, len( distr ) +  1 ), ybest, abestb, abestp )

psnum = -1

with open( "hdim_sph_s_table.txt" ) as f:
	for line in f:
		ls = line.split( "\t" )
		if ls[0] == "1":
			slist = ls[ 2: ]
			continue
		snum = snkeys.index( ls[0] )
		if snum != psnum:
			if psnum != -1:
				if len( yval ) > 4:
					(al,bl,rl,ttl,stderrl) = stats.linregress( xmval, yval )
					#fit,fstats = np.polynomial.polynomial.polyfit( xmval, yval, 1, full = True )#, w = rcval )
					if ttl < 0.1:
						ax.plot( [ 0, 9 ], [ bl, bl + 9 * al ], c="burlywood",lw=2 )
				cmax = max( cval )
				cs = ax.scatter( xval, yval, s = sval, c = cval, cmap = colmap, norm = mcolors.Normalize( vmin = 0, vmax = mmax ) )
			ax = plt.subplot( gs[ ( snum // 4, snum % 4 ) ] )
			ax.set_title( snames[ ls[0] ], loc = "left", size = "22" )
			if False:
				if ls[0] in [ "1118847", "840888" ]:
					ax.set_ylim( [ 0., 1. ] )
					ax.set_yticks( [ 0., 0.4, 0.8 ] )
				elif ls[0] in [ "72607", "848608" ]:
					ax.set_ylim( [ 0., 0.5 ] )
					ax.set_yticks( [ 0., 0.2, 0.4 ] )
				elif ls[0] in [ "217320" ]:
					ax.set_ylim( [ 0., 0.25 ] )
					ax.set_yticks( [ 0., 0.1, 0.2 ] )
				else:
					ax.set_ylim( [ 0., 0.125 ] )
					ax.set_yticks( [ 0., 0.1 ] )
			ax.set_xlim( [ -1, 10 ] )
			ax.set_xticks( [ 0, 1, 4, 5, 8, 9 ] )
			#ax.set_xticklabels( [ "L", "", "OV", "", "T", "" ] )
			ax.set_xticklabels( [ "d", "h", "d", "h", "d", "h",  ] )
			ax.tick_params( axis='y', labelsize = 18. )
			ax.tick_params( axis='x', labelsize = 18. )
			for ploc in range( 3 ):
				ax.text( ploc * 4 + 0.5, -0.08, [ "L", "OV", "T" ][ ploc ], va = "top", ha = "center", size = 22 ) 
				#ax.text( ploc * 4 + 1, -0.09, "h", va = "top", ha = "center", size = 18 ) 
			xval = []
			yval = []
			xmval = []
			cval = []
			rcval = []
			sval = []
			psnum = snum
		ax.set_ylim( [ 0, 1.5 ] )
		ax.set_yticks( [ 0., 0.5, 1 ] )
		dimres = ls[ 3 ]
		dval = dimres.split( "/" )
		if dval[2][0] == "-" or float( dval[2] ) < 0.8 or int( dval[5] ) == 0:
			continue
		distr = map( float, ls[4].strip().split() )
		( vx, vy, abest, abestl ) = calc_jakovenko( distr )		
		y = abestl
		percent = float( dval[4] ) / int( dval[5] )
		normed_abest = abest * ( mmax / pscales[ snum ] ) * ( len( distr ) / ( 0.05 - 0.003 ) ) * 0.8
		residuals = mmax - normed_abest
		#mmax = max( mmax, residuals )
		#mmin = max( mmin, residuals )
		sample = ls[2]
		loc = llet.index( sample[0] )
		hdis = 1 if sample in shealthy else 0
		x = 4 * loc + hdis
		xval.append( x )
		yval.append( y )
		cval.append( residuals ) #1. - percent )
		rcval.append( percent )
		xmval.append( 4 * loc )
		sval.append( 250 * math.sqrt( percent ) )

if len( yval ) > 4:
	(al,bl,rl,ttl,stderrl) = stats.linregress( xmval, yval )
	#fit,fstats = np.polynomial.polynomial.polyfit( xmval, yval, 1, full = True )#, w = rcval )
	if ttl < 0.1:
		ax.plot( [ 0, 9 ], [ bl, bl + 9 * al ], c="burlywood",lw=2 )
cmax = max( cval )
cs = ax.scatter( xval, yval, s = sval, c = cval, cmap = colmap, norm = mcolors.Normalize( vmin = 0, vmax = mmax ) )
		
if True:
	cbaxis = plt.subplot( gs[ (1,3) ] )
	cbar = fig.colorbar( cs, cax = cbaxis, orientation = "vertical", ticks = [ mmax, mmax * 0.6, mmax * 0.2 ] )
	psizes = [ 1, 5, 10, 20, 50 ]
	ylim = cbar.ax.get_ylim()
	cbar.ax.scatter( [ -1.5 ] * len( psizes ), [ ylim[0] + ( r + 1 ) * ( ylim[1] -  ylim[0] ) / ( len( psizes ) + 1 ) for r in range( len( psizes ) ) ], s = [ 25 * math.sqrt( psize ) for psize in psizes ], c = "#b9754b" )
	for r in range( len( psizes ) ):
		cbar.ax.text( -2,  ylim[0] + ( r + 1 ) * ( ylim[1] -  ylim[0] ) / ( len( psizes ) + 1 ), str( psizes[r] ) + " %", ha = "right", va = "center", size = 20 )
		
		
	#cbar.ax.set_yticks( [ mmax, mmax * 0.6, mmax * 0.2  ] )
	cbar.ax.set_yticklabels( [ "0", "b0 / 2", "b0"  ] )
	cbar.ax.tick_params( axis='y', labelsize = 20. )
	cbar.ax.set_xlim( [ -4, 1.2 ] )
#	cbar.ax.set_title( "Boltzmann rate", size = 20., loc="right" )
	
plt.savefig( "fig6.png" )

			
		
		

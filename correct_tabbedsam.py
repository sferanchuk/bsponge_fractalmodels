
import sys

with open( sys.argv[1] ) as f:
	for line in f:
		ls = line.strip().split( "\t" )
		lsm = ls[0] + "\t" + "\t".join( ls[-18:] )
		print lsm
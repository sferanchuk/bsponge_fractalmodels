while read line
do
	name=${line%.fa}
	~/soft/usearch10.0.240_i86linux32 -usearch_global ${name}.fa -db fdimrefs.udb -id 0.95 -samout gf_${name}.sam -strand both -maxaccepts 0 -maxrejects 0 -top_hit_only
	python2 correct_tabbedsam.py gf_${name}.sam >tmp.sam
	mv tmp.sam gf_${name}.sam
done <hdimsamples.lst


#iterate over the outputs from normalize-hic-segments.r
 for f in *.txt; do java -jar ../juicer/juicer_tools/juicer_tools_1.14.08.jar pre $f ${f%.txt}.hic hg38 -r 500000 -d;done
 
 chrs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
 
 for f in *.hic; do
	for chr in ${chrs[@]};do
		java -jar ../juicer/juicer_tools/juicer_tools_1.14.08.jar eigenvector KR $f $chr BP 500000 ../pca/${f%.hic}-pc1-$chr.txt
	done
 done
 
 
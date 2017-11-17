#! /bin/bash
"./hw3" > imm.txt \
	&& cat imm.txt | grep -o -E "1[0-9]+x[ATCG-]+x[ATCG-]+" | cut -d '1' -f2 | tr 'x' '\n' > naiveGap.txt \
	&& cat imm.txt | grep -o -E "2[0-9]+x[ATCG-]+x[ATCG-]+" | cut -d '2' -f2 | tr 'x' '\n' > affineGap.txt \
	&& rm imm.txt \
	&& rm seq.txt

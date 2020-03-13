#!/bin/zsh

md=( "" "." )
ed=( "" "." )

zparseopts -D -K -prefix:=pf -mapdir:=md -errdir:=ed

base=${pf[2]}

echo -n "#sample\tall_reads\ttrimmed_reads\tdup_rate\tmapped_reads"
echo -n "\tassigned_to_fragment\tpaired_end"
echo "\ttrue_links\tsame_fragment\tbetween_scaffolds\twithin_scaffolds"



echo -n $base"\t"
cat ${ed[2]}/${base}_cutadapt.err | grep -m 1 "^Total read pairs processed:" \
 | awk -F: '{printf $2"\t"}' | tr -d ', '
cat ${ed[2]}/${base}_cutadapt.err | grep -m1 'Pairs written (passing filters):' \
 | awk -F: '{print $2}' | tr -d ', ' | cut -d '(' -f 1 | tr '\n' '\t'
grep -m1 Improper ${ed[2]}/${base}_novosort1.err | tr -s ' ' | cut -d ' ' -f 4,5 | tr -d , \
 | awk '{printf $2/($2+$1)"\t"}' 
awk '{printf $1"\t"}' ${md[2]}/${base}_both_mapped_q10.len
awk '{printf $1"\t"}' ${md[2]}/${base}_reads_to_fragments.bed.len 
awk '{printf $1"\t"}' ${md[2]}/${base}_pe_count.txt
cat ${md[2]}/${base}_frag_stat.txt


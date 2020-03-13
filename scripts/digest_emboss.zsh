#!/bin/zsh

# Perform an in silico digest of an assembly FASTA

e=(-restrict '/opt/Bio/EMBOSS/6.6.0/bin/restrict')
t=(-bedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools')
b=(-rebase '/filer-dg/agruppen/DG/mascher/rebase_181023/link_emboss_e.txt')

zparseopts -D -K -- -bedtools:=t -restrict:=e -rebase:=b -name:=n -sitelen:=l -enzyme:=s -ref:=r -minlen:=m 

bedtools=$t[2]
ref=$r[2]
name=$n[2] #not actually used ... ?
emboss=$e[2]
enzyme=$s[2]
minlen=$m[2]
len=$l[2]
rebase=$b[2]

bed=${ref:r}_${enzyme}_fragments_${minlen}bp.bed

{
 $emboss -enzymes $enzyme -sitelen $len -sequence $ref -datafile $rebase -outfile /dev/stdout | awk NF \
  | awk '/^# Sequence:/ {printf "\n"$3} /^ *[0-9]/ {printf "\t"0+$1}' | awk NF \
  | awk -v len=$len '{for(i = 2; i < NF; i++) print $1"\t"$i+len-1"\t"$(i+1)-1}' \
  | sort -S10G -k 1,1 -k 2,2n \
  | awk '$3 - $2 >= '$minlen > $bed && \
 awk 'BEGIN{OFS="\t"}
      $3 - $2 <= 400 {print $0, $1":"$2"-"$3; next} 
      {print $1, $2, $2+200, $1":"$2"-"$3;
       print $1, $3-200, $3, $1":"$2"-"$3}' $bed \
  | $bedtools nuc -bed - -fi $ref > ${bed:r}_split.nuc.txt 
} 2> ${ref:r}_${enzyme}_fragments_${minlen}bp_digest.err  

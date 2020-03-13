#!/bin/zsh

add_dist='
BEGIN{
 OFS="\t"
}

$9 == "-" {
 d1=$2-$12
}

$10 == "-" {
 d2=$5-$15
}

$9 == "+" {
 d1=$13-$3
}

$10 == "+" {
 d2=$16-$6
}

d1 < 0 {
 d1 = 0
}

d2 < 0 {
 d1 = 0
}

$1 == $4 {
 if($12 == $15)
  rd = 0
 else if($13 < $15)
  rd = $15-2-($13+3)
 else
  rd = $12-2-($16+3)
}

$1 != $4 {
 rd=-1
}

{
 l1=d1+$3-$2
 l2=d2+$6-$5
 print $0"\t"d1,d2,l1,l2,l1+l2,rd
}'

dist_stat='
BEGIN{
 OFS="\t"
}

#same fragment
$22 == 0 {
 same++
}

#different seq scaffold
$22 == -1 {
 diffcl++
}

#adj frags
#$22 == 1 {
# adj++
#}

#same seq scaffold
$22 > 1 {
 samecl++
}

END {
 print NR,0+same,0+diffcl,0+samecl
}'

zparseopts -D -K -- -linker:=l -data_dir:=d -err_dir:=e -out_dir:=O -tmp:=p -name:=n \
 -ref:=r -bed:=b -minlen_read:=m -maxlen:=x -threads:=T -minq:=q \
 -batchmem:=bm -indexsize:=is \
 -onlyfrag:=o \
 -cutadapt:=cut -minimap:=mini -samtools:=sam \
 -bedtools:=bedt -novosort:=novo -bgzip:=bgz

bmem=$bm[2]
isize=$is[2]
out_dir=$O[2]
err_dir=$e[2]
bgzip=$bgz[2]
novosort=$novo[2]
bedtools=$bedt[2]
samtools=$sam[2]
minimap=$mini[2]
cutadapt=$cut[2]
linker=$l[2]
minlen=$m[2]
threads=$T[2]
ref=$r[2]
tmp=$p[2]
minq=$q[2]
maxlen=$x[2]
bed=$b[2]
onlyfrag=$o[2]
data_dir=$d[2]
name=$n[2]


if [[ $threads -gt 16 ]]; then
 novosort_threads=16
else
 novosort_threads=$threads
fi

name=$name
bam=${name}.bam
minimaperr=$err_dir/${name}_minimap.err
samtoolserr=$err_dir/${name}_samtools.err
sorterr1=$err_dir/${name}_novosort1.err
sorterr2=$err_dir/${name}_novosort2.err
cutadapterr=$err_dir/${name}_cutadapt.err

b=$name:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

# Note the find commands. If sample names contain R1 and R2, this will not work.
if [[ $onlyfrag -eq 0 ]]; then
 $cutadapt -f fastq --interleaved -a $linker -A $linker -O 1 -m $minlen 2> $cutadapterr \
   <(find $data_dir/*_R1*.f*.gz | sort | xargs zcat) \
   <(find $data_dir/*_R2*.f*.gz | sort | xargs zcat) \
 | $minimap -ax sr -R $rgentry -t $threads -2 $ref /dev/stdin 2> $minimaperr \
 | $samtools view -Su /dev/stdin 2> $samtoolserr \
 | $novosort -c $novosort_threads -t $tmp -m $bmem --keepTags --md /dev/stdin 2> $sorterr1 \
 | $novosort -c $novosort_threads  -t $tmp -m $bmem -n -o $bam /dev/stdin 2> $sorterr2
 
echo $pipestatus > $err_dir/${name}_pipestatus.txt
fi


f="${name}_reads_to_fragments.bed.gz"
err="${name}_reads_to_fragments.err"

{
 {
  $samtools view -H $bam
  $samtools view -q 10 -F 3332 $bam \
   | tee >(cut -f 1 | uniq -c | awk '$1 == 2' | wc -l > ${name}_both_mapped_q10.len) \
   | cut -f -16 \
   | awk 'old != $1 {old=$1; printf "\n"$0; next} {printf "\t"$0}' \
   | awk -F '\t' 'NF == 32' \
   | awk '{for(i=1; i<=15; i++) printf $i"\t"; printf $16"\n"$17; for(i=18; i<=32; i++) printf "\t"$i; print ""}'
 } | $samtools view -Su - \
  | $bedtools pairtobed -f 1 -type both -abam - -bedpe -b $bed \
  | awk 'old != $7 {old=$7; printf "\n"$0; next} {printf "\t"$0}' \
  | awk 'NF == 26'  | cut -f 1-13,24-26 | awk "$add_dist" \
  | tee >($bgzip -c > $f) \
        >(wc -l > $f:r.len) \
	>(awk -v maxlen=$maxlen '$21 > maxlen && $1 == $4' \
	   | cut -f 1,2,3,5,6 | awk '{print $5 - $2}' \
	   | awk -v maxlen=$maxlen '$1 <= maxlen' | sort | uniq -c > ${name}_length_dist_PE.txt ) \
        >(awk -v maxlen=$maxlen '$21 > maxlen' | wc -l > ${name}_pe_count.txt) \
	>(awk -v maxlen=$maxlen '$21 <= maxlen' | cut -f 21 | sort | uniq -c > ${name}_length_dist.txt) \
        >(awk -v maxlen=$maxlen '$21 <= maxlen' | awk "$dist_stat" > ${name}_frag_stat.txt) \
 | awk -v maxlen=$maxlen '$21 <= maxlen' \
 | cut -f 11,12,14,15 | awk '$1 != $3 || $2 != $4' \
 | awk '{print $0"\n"$3"\t"$4"\t"$1"\t"$2}' \
 | $bgzip > ${name}_fragment_pairs.tsv.gz
} 2> $err

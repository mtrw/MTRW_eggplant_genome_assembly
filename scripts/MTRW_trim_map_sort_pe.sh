#!/usr/bin/bash 


seqtk='/opt/Bio/seqtk/1.0-r76/bin/seqtk'
bbduk='/opt/Bio/bbmap/37.28/bin/bbduk.sh'
minimap2='/home/wallace/bin/minimap2/minimap2'
novosort='/opt/Bio/novocraft/V3.06.05/bin/novosort'
samtools='/opt/Bio/samtools/1.9/bin/samtools'

adapt1='AGATCGGAAGAGC'
adapt2='AGATCGGAAGAGC'
indexsize='50G'
mem='50G'
threads=10
minimap_preset='sr'

while getopts ':r:q:a:A:p:' opt $@; do
	case ${opt} in
		r )
			ref_idx=$OPTARG
		;;
		q )
			reads_dir=$OPTARG
		;;
		a )
			adapt1=$OPTARG
		;;
		A )	
			adapt2=$OPTARG
		;;
		p )
			minimap_preset=$OPTARG
		;;
	esac;
done

base=${reads_dir##*/}
rgentry="@RG\tID:${base}\tPL:ILLUMINA\tPU:${base}\tSM:$base"
out_bam=${base}.bam
temp_dir=${RANDOM}_tmp
mkdir ${temp_dir}

printf "	Ref: $ref_idx
	Reads in: $reads_dir
	Adapters: $adapt1 $adapt2
	Minimap setting: $minimap_preset
	Regentry: ${rgentry}
	Temp dir: ${temp_dir}\n" 1>&2


if ${samtools} quickcheck $out_bam; then
	(>&2 echo $out_bam is already in tip top shape: exiting ...)
	exit 0
else
	(>&2 echo Beginning mapping for $out_bam)
fi

${seqtk} mergepe  <( zcat $( find  ${reads_dir}/*_R1_*.f* ) ) <( zcat $( find  ${reads_dir}/*_R2_*.f* ) ) | \
${bbduk} \
	in1='stdin.sh' \
	interleaved='t' \
	literal=${adapt1},${adapt2} \
	out='stdout.fq' \
	ktrim='r' \
	k=${#adapt1} \
	mink=2 \
	qtrim='rl' \
	trimq=20 \
		2> ${base}_bbduk.err | \
${minimap2} -a -x ${minimap_preset} -R ${rgentry} -t $threads -2 -I $indexsize -K $mem $ref_idx /dev/stdin \
	2> ${base}_minimap.err | \
${samtools} view -b /dev/stdin \
	2> ${base}_samtools.err | \
${novosort} -c $threads -t ${temp_dir} -m $mem -i --keepTags --md -o $out_bam /dev/stdin \
	2> ${base}_novosort.err


echo $pipestatus | tr ' ' '\n' | grep -q '^[^0$]' || (>&2 echo "$base: No errors detected ... ")

rm -r ${temp_dir}

exit 0



#!/bin/awk -f

BEGIN{
 OFS="\t"
 h["1/1"]=2
 h["0/1"]=1
 h["0/0"]=0
}

/^#CHROM/ && !header {
 for(i = 1; i <= NF; i++)
  s[i]=$i
 header=1
}

/INDEL/ || $4 == "N" || $5 ~ /,/ || /^#/ { 
 next
}


{

 printf $2"---"$2"---"$4"---"$5"---"$6"---"s[i]"---"h[a[1]]"---"a[5]"---"a[3]"---"a[4]"\n"
 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  #if(a[3] > 0)
   print $1, $2, $4, $5, $6, s[i], h[a[1]], a[5], a[3], a[4]
 }
}

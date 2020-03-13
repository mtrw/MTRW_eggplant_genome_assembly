library(igraph)
library(zoo)

make_hic_info<-function(cluster_info, super_global, chrs){
 s<-super_global$super_info
 s[!duplicated(s$chr),]->s
 s[chr %in% chrs]->s

 super_global$membership[, .(cluster, super, bin, rank, backbone)]->tmp
 tmp[super %in% s$super]->tmp
 tmp[, super := NULL]
 setnames(tmp, c("cluster", "hic_bin", "hic_rank", "hic_backbone"))
 tmp[cluster_info, on="cluster"]->cluster_info
 cluster_info[order(chr, hic_bin, hic_rank, cluster)]
}

inspector_export<-function(hic_map, assembly, inversions, species, file){
  hic_map$hic_1Mb$norm->x
  wheatchr(agp=T, species=species)[, .(chr1=chr, agp_chr)][x, on="chr1"]->x
  x[, c("chr1", "chr2", "id2", "id1", "dist") := list(NULL, NULL, NULL, NULL, NULL)]
  
  inversions$ratio[, .(agp_chr, bin, r)] -> y
  
  copy(hic_map$agp)[, chr := NULL][]->a
  assembly$info[, .(scaffold, popseq_chr=popseq_alphachr)][a, on="scaffold"]->a
  saveRDS(list(agp=a, binsize=1e6, links=x, ratio=y), file=file)
}

insert_node<-function(df, el){
 while({
  setkey(el[cluster2 %in% df$cluster & cluster1 %in% unique(el[!cluster1 %in% df$cluster & cluster2 %in% df$cluster]$cluster1)], "cluster2")[setkey(setnames(df[, .(cluster, bin)], c("cluster2", "bin")), "cluster2"), allow.cartesian=T][!is.na(cluster1)][
   order(cluster1, bin)][, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]->y
   y[order(cluster1, bin)]->y
   length(which(y$dist == 1)->idx) > 0
  }) {

   setkeyv(setnames(el[, .(cluster1, cluster2, weight)], c("path_node1", "path_node2", "old_path")), c("path_node1", "path_node2"))[
    setkeyv(data.table(cluster=y[idx, cluster1], path_node1=y[idx, cluster2], path_node2=y[idx+1, cluster2], 
	      weight1=y[idx, weight], weight2=y[idx+1, weight], bin=y[idx, bin]), c("path_node1", "path_node2"))]->z
    head(z[ ,diff:= weight1 + weight2 - old_path][order(diff)],1)->z
   
   m=z$bin
   data.table(rbind(df[1:m], data.table(cluster=z$cluster, bin=m+1), df[(m+1):nrow(df),][, bin:=bin+1][]))->df
 }
 df
}

node_relocation<-function(df, el, maxiter=100, verbose=T){
 i<-0
 while({
  i<-i+1
  df[order(bin)]->df
  x<-data.table(old_node1=df[1:(nrow(df)-2)]$cluster, cluster=df[2:(nrow(df)-1)]$cluster, old_node2=df[3:nrow(df)]$cluster)
  setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
  setnames(setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))[setkeyv(x, c("old_node1", "old_node2"))], "weight", "new_edge3")->x
  setnames(setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))[setkeyv(x, c("cluster", "old_node1"))], "weight", "old_edge1")->x
  setnames(setnames(ee, "old_node1", "old_node2")[setkeyv(x, c("cluster", "old_node2"))], "weight", "old_edge2")->x
  x[!is.na(new_edge3) ]->x

  setkey(el[, .(cluster1, cluster2)], "cluster2")[setkey(copy(df), "cluster")][order(cluster1, bin)]->t
  which(t[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]$dist == 1)->idx

  data.table(cluster=t$cluster1[idx], new_node1=t$cluster2[idx], new_node2=t$cluster2[idx+1])->t
  setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
  setnames(setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))[setkeyv(t, c("cluster", "new_node1"))], "weight", "new_edge1")->t
  setnames(setnames(ee, "new_node1", "new_node2")[setkeyv(t, c("cluster", "new_node2"))], "weight", "new_edge2")->t
  setnames(setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))[setkeyv(t, c("new_node1", "new_node2"))], "weight", "old_edge3")->t
  setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]->m
  m[ ,diff := new_edge1+new_edge2+new_edge3 - old_edge1 - old_edge2 - old_edge3]
  nrow(m<-m[diff < 0][order(diff)]) > 0 & i <= maxiter
  }) {
  ne<-head(data.frame(m), 1)
  idx<-data.frame(cluster=df$cluster, idx=1:nrow(df))
  invisible(lapply(c("cluster", "new_node1", "new_node2", "old_node1", "old_node2"), function(i) {
   merge(ne, by.x=i, by.y="cluster", idx)->>ne
   colnames(ne)[which(colnames(ne) == "idx")]<<-paste(sep="_", "idx", i)
  }))
  df[with(ne, {
   min_new <- min(idx_new_node1, idx_new_node2)
   max_new <- max(idx_new_node1, idx_new_node2)
   min_old <- min(idx_old_node1, idx_old_node2)
   max_old <- max(idx_old_node1, idx_old_node2)
   if(min_old < min_new) {
    c(1:min_old, max_old:min_new, idx_cluster, max_new:nrow(df))
   } else {
    c(1:min_new, idx_cluster, max_new:min_old, max_old:nrow(df))
   }
  }),]->df
  df$bin<-1:nrow(df)
  df<-data.table(df)
 }
 if(verbose){
  cat(paste0("Node relocation steps: ", i-1, "\n"))
 }
 df
}

kopt2<-function(df, el){
 setkeyv(el, c("cluster1", "cluster2"))->el
 while({
  df[, .(cluster1=cluster[1:(nrow(df)-1)], cluster2=cluster[2:nrow(df)])]->m
  el[setkeyv(m, c("cluster1", "cluster2"))]->m
  m[ ,.(cluster1, cluster2, weight)]->m
  setnames(m, "weight", "weight12")->m
  m[, dummy:=1]
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster1", "bin1")), "cluster1")[setkey(m, "cluster1")]->m
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster2", "bin2")), "cluster2")[setkey(m, "cluster2")]->m
  copy(m)->n
  setnames(n, c("cluster1", "cluster2", "bin1", "bin2", "weight12"), c("cluster3", "cluster4", "bin3", "bin4", "weight34"))
  setkey(m, "dummy")[setkey(n, "dummy"), allow.cartesian=T]->mn
  mn[, dummy:=NULL]
  mn[bin1 < bin3]->mn
  o<-el[, .(cluster1, cluster2, weight)]
  setkeyv(setnames(copy(o), c("cluster1", "cluster3", "weight13")), c("cluster1", "cluster3"))[setkeyv(mn, c("cluster1", "cluster3"))]->mn
  setkeyv(setnames(copy(o), c("cluster2", "cluster4", "weight24")), c("cluster2", "cluster4"))[setkeyv(mn, c("cluster2", "cluster4"))]->mn
  mn[, old:=weight12+weight34]
  mn[, new:=weight13+weight24]
  mn[, diff:=old-new]
  mn[order(-diff)]->mn
  nrow(mn[diff > 0]) > 0 
 }) {
  head(mn, 1)->x
  bin1<-df[cluster == x$cluster1]$bin
  bin2<-df[cluster == x$cluster2]$bin
  bin3<-df[cluster == x$cluster3]$bin
  bin4<-df[cluster == x$cluster4]$bin
  df[c(1:bin1, bin3:bin2, bin4:nrow(df))]->df
  df[, bin:=1:nrow(df)]->df
 }
 df
}

make_super<-function(hl, cluster_info, prefix="super", cores=1, paths=T, path_max=0, known_ends=F, 
		     maxiter=100, verbose=T){
 cat("In make super routine \n")
 hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl

 hl[cluster1 < cluster2]->e
 graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
 E(g)$weight<-e$weight

 data.table(cluster=V(g)$name, super=paste(prefix, sep="_", clusters(g)$membership))->mem
 cluster_info[mem, on="cluster"]->mem

 mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))), keyby=super]->info
 mem[, .(cluster1=cluster, super)][hl, on="cluster1"]->e

 list(super_info=info, membership=mem, graph=g, edges=e)->s

 if(paths){
  if(path_max > 0){
   idx<-head(s$super_info[order(-length)], n=path_max)$super
  } else {
   idx<-s$super_info$super
  }
  rbindlist(mclapply(mc.cores=cores, idx, function(i) {
   start <- end <- NULL
   if(known_ends){
    s$mem[super == i & !is.na(cM)][order(cM)]$cluster->x
    start=x[1]
    end=tail(x,1)
   }
   make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)->x
   if(verbose){
    cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
   }
   x
  }))[s$membership, on="cluster"]->s$membership
 }

 s
}

make_super_path<-function(super, idx=NULL, start=NULL, end=NULL, maxiter=100, verbose=T){
 write("In make super path routine",stderr())
 submem<-super$mem[super == idx]
 super$edges[super == idx, .(cluster1, cluster2, weight)]->el
 #if(submem$agp_chr[1]==9){browser()}
 minimum.spanning.tree(induced.subgraph(super$graph, submem$cluster))->mst
 E(mst)$weight <- 1

 if(is.null(start) | is.null(end)){
  V(mst)[get.diameter(mst)]$name->dia
 } else {
  V(mst)[get.shortest.paths(mst, from=start, to=end)$vpath[[1]]]$name->dia
 }

 data.table(cluster=dia, bin=1:length(dia))->df

 df<-insert_node(df, el)
 df<-kopt2(df, el)
 df<-node_relocation(df, el, maxiter=maxiter, verbose=verbose)

 data.frame(df)->df
 data.frame(el)->el
 data.frame(cluster=df$cluster, rank = 0)->ranks
 r=0

 while(length(n<-unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1)) > 0) {
  r = r+1
  subset(el, cluster2 %in% df$cluster & cluster1 %in% n)->tmp
  tmp[!duplicated(tmp$cluster1),]->tmp
  rbind(ranks, data.frame(cluster=tmp$cluster1, rank=r))->ranks
  merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
  rbind(df, data.frame(cluster=x$cluster1, bin=x$bin))->df
 }

 merge(df, submem)->df
 df$bin<-as.integer(df$bin)
 flip<-with(df, suppressWarnings(cor(bin, cM))) < 0
 if((!is.na(flip)) & flip) {
  with(df, max(bin) - bin + 1)->df$bin
 }
 ranks$rank<-as.numeric(ranks$rank)
 merge(df, ranks)->df
 df$backbone <- df$cluster %in% dia
 data.table(df)[, .(cluster, bin, rank, backbone)]
}

make_hic_map<-function(hic_info, links, ncores=1){
 cat("In make hic map routine \n")
 copy(links)->hl
 copy(hic_info)->info
 
 setnames(hl, c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))
 setnames(info, "scaffold", "cluster")
 setkey(info, "cluster")
 info[!is.na(chr), length(unique(chr))]->nchr
 
 make_hic_info(info, 
               
 super_global<-make_super(hl, cluster_info=info, cores=ncores, known_ends=T, path_max=nchr), chrs=1:nchr)->res
 res[order(chr, hic_bin)][, .(scaffold=cluster, chr, cM, hic_bin, hic_backbone, hic_rank)][!is.na(hic_bin)]
}

correct_inversions<-function(hic_map, scaffolds, species="wheat"){
 cat("In correct inversions routine, inverting:\n")
 copy(hic_map$hic_map)->z
 print(z[scaffold %in% scaffolds]$scaffold)
 z[scaffold %in% scaffolds, consensus_orientation := ifelse(is.na(consensus_orientation), -1, -1 * consensus_orientation)]
 make_agp(z, gap_size=hic_map$gap_size, species=species)->a
 
 list()->new
 new$hic_map <- z
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_inversions <- scaffolds
 new
}

hic_map<-function(info, assembly, frags, species="wheat", ncores=1, min_scaff_length=0 , min_nfrag_scaffold=50, max_cM_dist = 20, 
		  binsize=5e5, min_nfrag_bin=30, gap_size=100, orient=T, agp_only=F, map=NULL){
  cat("In hic map routine \n")
 if(!agp_only){
  copy(info)->hic_info
  short_scaffs <- assembly$info[length < min_scaff_length]$scaffold
  hic_info[, excluded:=(nfrag < min_nfrag_scaffold | scaffold %in% short_scaffs) ] #per scaffold, number of restriction frags, popseq position, chromosome assignment.
  assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl #interscaffold fragment pair hic link counts (symmetrical: every entry scaff1 --> scaff2 is accompanied by an entry scaff2 --> scaff1)
  hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl #attach info to each inter-scaffold link count (row) for each scaffold in the link
  hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl #for the second scaffold each link
  hl[chr1 == chr2]->hl #only intrachromosomal links #this excludes NA chromosomes and gets rid of our broken scaffolds! ######################
  hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)] #only close (popseq-wise) links (or those with no popseq info)
  hl[, weight:=-log10(nlinks)]

  cat("Scaffold map construction started.\n")
  make_hic_map(hic_info=hic_info, links=hl, ncores=ncores)->hic_map
  cat("Scaffold map construction finished.\n")

  if(orient){
   options(scipen = 1000)
   frags[, .(n=.N), keyby=.(scaffold, pos = start %/% binsize * binsize)]->fragbin
   fragbin[, id := paste(sep=":", scaffold, pos)]
   fragbin<-hic_info[excluded == F][fragbin, on="scaffold", nomatch=0]
   fragbin[, excluded := NULL]

   assembly$fpairs[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize * binsize, scaffold2, pos2 = pos2 %/% binsize * binsize)]->binl
   binl[, id1 := paste(sep=":", scaffold1, pos1)]
   binl[, id2 := paste(sep=":", scaffold2, pos2)]
   binl[id1 != id2]->binl
   fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
   fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
   binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
   setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))
   
   cat("Scaffold bin map construction started.\n")
   fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
   hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
   binl[chr1 == chr2 & (abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2))]->binl
   binl[, weight:=-log10(nlinks)]
   make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores)->hic_map_bin
   cat("Scaffold bin map construction finished.\n")
   w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
   w<-w[, .(gbin=mean(na.omit(hic_bin)),
	    hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
   hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0

   short_lgs <- z0[,.N,by="chr"][N<3]$chr
   z0 <- z0[!chr %in% short_lgs] #ADDED TO MAKE IT SKIP 2-SCAFFOLD LINKAGE GROUPS ("CHROMOSOMES")
   z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
   z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
   z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
   z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
   w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
   w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
   w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
   z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
		   suppressWarnings(cor(x[1:3], x[4:6]))
		   })]
   z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
   ccor[w]->m
   m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
   m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented

   setnames(hic_map_oriented, "chr", "consensus_chr")
   setnames(hic_map_oriented, "cM", "consensus_cM")
   hic_map_oriented[, consensus_orientation := hic_orientation]
  } else {
   hic_map_oriented<-copy(hic_map)
   setnames(hic_map_oriented, "chr", "consensus_chr")
   setnames(hic_map_oriented, "cM", "consensus_cM")
   hic_map_oriented[, consensus_orientation := as.numeric(NA)]
   hic_map_oriented[, hic_cor := as.numeric(NA)]
   hic_map_oriented[, hic_invert := as.numeric(NA)]
   hic_map_oriented[, hic_orientation := as.numeric(NA)]
   hic_map_bin <- NA
  }

  if("orientation" %in% names(hic_info)){
   hic_info[, .(scaffold, old_orientation=orientation)][hic_map_oriented, on="scaffold"]->hic_map_oriented
   hic_map_oriented[!is.na(old_orientation), consensus_orientation := old_orientation]
   hic_map_oriented[, old_orientation := NULL]
  }

  hic_map_oriented[assembly$info, on="scaffold"]->hic_map_oriented
 } else {
  hic_map_oriented <- map$hic_map
  hic_map_bin <- map$hic_map_bin
  min_nfrag_scaffold <- map$min_nfrag_scaffold
  binsize <- map$binsize
  max_cM_dist <- map$max_cM_dist
  min_nfrag_bin <- map$min_nfrag_bin
  gap_size <- map$gap_size
 }

 make_agp(hic_map_oriented, gap_size=gap_size, species=species)->a

 a$agp[, .(length=sum(scaffold_length)), key=agp_chr]->chrlen
 chrlen[, alphachr := sub("chr", "", agp_chr)]
 chrNames(species=species)[chrlen, on="alphachr"]->chrlen
 chrlen[, truechr := !grepl("Un", alphachr)]
 chrlen[order(!truechr, chr)]->chrlen
 chrlen[, offset := cumsum(c(0, length[1:(.N-1)]))]
 chrlen[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]

 list(agp=a$agp, agp_bed=a$agp_bed, chrlen=chrlen, hic_map=hic_map_oriented, hic_map_bin=hic_map_bin)->res
 invisible(lapply(sort(c("min_nfrag_scaffold", "max_cM_dist", "binsize", "min_nfrag_bin", "gap_size")), function(i){
  res[[i]] <<- get(i)
 }))
 res
}

make_agp<-function(hic_map_oriented, gap_size=100, species="wheat"){
  cat("In make agp routine \n")
 #browser()
 hic_map_oriented[, .(scaffold, chr = consensus_chr,
	  popseq_cM=cM, # ifelse(consensus_chr == popseq_chr | is.na(consensus_chr), popseq_cM, NA),
	  scaffold_length = length, hic_bin, orientation=consensus_orientation)]->z
 chrNames(agp=T)[z, on="chr"]->z

 z[, agp_chr := "chrUn"]
 z[!is.na(hic_bin), agp_chr := sub("NA", "Un", paste0("chr", alphachr))]
 z[, alphachr := NULL]
 z[order(agp_chr, hic_bin, chr, popseq_cM, -scaffold_length)]->z
 z[, index := 2*1:.N-1]
 z[, gap := F]
 rbind(z, data.table(scaffold="gap", gap=T, chr=NA, popseq_cM=NA, scaffold_length = gap_size, hic_bin = NA, orientation = NA, agp_chr=z$agp_chr, index=z$index+1))->z
 z[order(index)][, head(.SD, .N-1), by=agp_chr]->z
 z[, agp_start := cumsum(c(0, scaffold_length[1:(.N-1)]))+1, by = agp_chr]
 z[, agp_end := cumsum(scaffold_length), by = agp_chr]

 z[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr=agp_chr)]->agp_bed

 list(agp=z, agp_bed=agp_bed)
}

read_fragdata<-function(info, map_10x=NULL, assembly_10x=NULL, file){
 fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
 fragbed[, length := end - start]
 fragbed[, start := start + 1]
 info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
 fragbed[, start := start - orig_start + 1]
 fragbed[, end := end - orig_start + 1]
 fragbed[, orig_start := NULL]
 fragbed[, orig_scaffold := NULL]
 if(!is.null(assembly_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][fragbed, on="scaffold"]->fragbed
  fragbed[orientation == 1, start := super_start - 1 + start]
  fragbed[orientation == 1, end := super_start - 1 + end]
  fragbed[orientation == -1, start := super_end - end + 1]
  fragbed[orientation == -1, end := super_end - start + 1]
  fragbed[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(fragbed, "super", "scaffold")
  fragbed[, .(nfrag = .N), keyby=scaffold][assembly_10x$info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 } else {
  fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 }
 list(bed=fragbed[], info=z[])
}

add_psmol_fpairs<-function(assembly, hic_map, nucfile, map_10x=NULL, assembly_10x=NULL, cov=NULL){
 if(is.null(map_10x)){
  assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 } else {
  assembly_10x$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 }
 hic_map$agp[agp_chr != "chrUn", .(chr, scaffold, orientation=orientation, agp_start, agp_end)]->a
 a[is.na(orientation), orientation := 1]
 setnames(copy(a), paste0(names(a), 1))[z, on="scaffold1"]->z
 setnames(copy(a), paste0(names(a), 2))[z, on="scaffold2"]->z
 z[orientation1 == 1, pos1 := agp_start1 - 1 + pos1]
 z[orientation1 == -1, pos1 := agp_end1 + 1 - pos1]
 z[orientation2 == 1, pos2 := agp_start2 - 1 + pos2]
 z[orientation2 == -1, pos2 := agp_end2 + 1 - pos2]
 z[!is.na(chr1) & !is.na(chr2), .(chr1, chr2, start1=pos1, start2=pos2)]->links

 n<-c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
 nuc<-fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
 assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
 z[, start := orig_start - off + 1]
 z[, end := orig_end - off + 1]
 if(!is.null(map_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
  z[orientation == 1, start := super_start - 1 + start]
  z[orientation == 1, end := super_start - 1 + end]
  z[orientation == -1, start := super_end - end + 1]
  z[orientation == -1, end := super_end - start + 1]
  z[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(z, "super", "scaffold")
 }
 a[z, on="scaffold", nomatch=0]->z
 z[orientation == 1, start := agp_start - 1 + start]
 z[orientation == -1, start := agp_end + 1 - start]
 z[orientation == 1, end := agp_start - 1 + end]
 z[orientation == -1, end := agp_end + 1 - end]
 z[, c("chr", "start", "end", n[4:length(n)]), with=F]->z
 z[, cov := 1]
 z->frags
 
 copy(hic_map) -> res
 res$links <- links
 res$frags <- frags
 res
}

contact_matrix<-function(hic_map, links, file, species="wheat", chrs=NULL, boundaries=T, grid=NULL, ncol=100, trafo=NULL, v=NULL){
 colorRampPalette(c("white", "red"))(ncol)->whitered

 if(is.null(chrs)){
  chrs <- unique(links$chr1) 
 }

 binsize <- min(links[dist > 0]$dist)

 links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
 if(is.null(trafo)){
  z[, col := whitered[cut(l, ncol, labels=F)]]
 } else {
  z[, col := whitered[cut(trafo(l), ncol, labels=F)]]
 }

 chrNames(agp=T, species=species)[z, on="chr"]->z

 pdf(file)
 lapply(chrs, function(i){
  z[chr == i, plot(0, las=1, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, main=chrNames(species=species)[chr == i, alphachr])]
  if(boundaries){
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
  }
  if(!is.null(grid)){
   max(hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr]]$agp_end)->end
   abline(v=seq(0, end, grid)/1e6, col="blue", lty=2)
  }
  if(!is.null(v)){
   abline(v=v/1e6, col="blue", lty=2)
  }
  z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
 })
 dev.off()
}

correct_5A<-function(hic_info, bam, map, assembly){
 fread(paste("samtools view -q 30 -F260", bam, "| cut -f 1,3,4"))->i90sam
 setnames(i90sam, c("i90_marker", "scaffold", "pos"))

 i90map<-fread(map, select=c(1,3,4))
 setnames(i90map, c("i90_marker", "i90_alphachr", "i90_cM"))
 i90map[, n := .N, by=i90_marker]
 i90map[n == 1][, n := NULL][]->i90map
 setnames(chrNames(species="wheat"), c("i90_alphachr", "i90_chr"))[i90map, on="i90_alphachr"]->i90map
 i90map[i90sam, on="i90_marker"][i90_chr == 5, .(orig_scaffold=scaffold, orig_pos=pos, cM=i90_cM)]->z
 copy(z)->i90k_5A

 assembly$info[, .(scaffold, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
 z[, .(scaffold, pos=orig_pos-orig_start+1, cM)]->i90k_5A

 z[, .(chr=5, icM=median(cM)), key="scaffold"][hic_info, on=c("scaffold", "chr")]->hic_info
 mx <- hic_info[chr == 5, max(na.omit(icM))]
 hic_info[chr == 5, cM := mx - icM]->hic_info
 hic_info[, icM := NULL][]
}

plot_inversions<-function(inversions, hic_map, bad=c(), species="wheat", file, height=5, width=20){
 yy<-inversions$summary
 w<-inversions$ratio
 yy[scaffold %in% bad]->yy2

 pdf(file, height=height, width=width)
 lapply(sort(unique(w$chr)), function(i){
  w[chr == i, plot(bin/1e6, r, pch=20, las=1, bty='l', ylab="r", xlab="genomic position (Mb)", main=alphachr[i])]
  hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i]$agp_chr, abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
  if(nrow(yy2) > 0 & i %in% yy2$chr){
   yy2[chr == i, rect(col="#FF000011", agp_start/1e6, -1000, agp_end/1e6, 1000), by=scaffold]
  }
 })
 dev.off()
}

find_inversions<-function(hic_map, links, species="wheat", chrs=NULL, cores=1, winsize=15, maxdist=1e8, threshold=40, factor=100){
  cat("In find inversions routine \n")
 if(is.null(chrs)){
  unique(links$chr1) -> chrs
 }

 links[chr1 %in% chrs]->zz
 zz[, l := log(factor*nlinks_norm/sum(nlinks_norm)*1e6)]->zz
 hic_map$chrlen[, .(chr1=chr, length)][zz, on="chr1"]->zz
 zz[dist <= bin1 & dist <= bin2 & (length - bin1) >= dist & (length - bin2) >= dist]->zz
 zz[l >= 0 & dist <= maxdist, .(r = sum(sign(bin1 - bin2) * l)), key=.(chr=chr1, bin=bin1)]->w
 chrNames(agp=T, species=species)[w, on="chr"]->w

 data.table(w, rbindlist(mclapply(mc.cores=cores, chrs, function(i){
  ww<-w[chr == i][order(bin)]
  data.table(
   cc=ww[, rollapply(r, width=winsize, FUN=function(x) cor(x, 1:length(x)), align="left", fill=NA)],
   lm=ww[, rollapply(r, width=winsize, FUN=function(x) lm(data=.(a=x, b=1:winsize/winsize), a~b)$coefficient[2], align="left",fill=NA)]
  )
 })))->w

 hic_map$agp[gap == F, .(agp_chr, bin=agp_start, scaffold)]->x
 copy(w)->y
 y[, bin := bin + 1]
 x[y, on=c("agp_chr", "bin"), roll=T]->y
 y[agp_chr != "chrUn"]->y
 y[, .(s=sd(r)), key=scaffold][!is.na(s)][order(-s)]->yy
 hic_map$hic_map[, .(scaffold, chr=consensus_chr, hic_bin)][yy, on="scaffold"]->yy
 hic_map$agp[, .(scaffold, agp_start, agp_end)][yy, on="scaffold"]->yy
 list(ratio=w, summary=yy)
}
 
correct_multi_inversions<-function(hic_map, ranges, species="wheat"){
 copy(ranges)->y
 chrNames(agp=T, species=species)[, .(consensus_chr=chr, agp_chr)][hic_map$hic_map, on="consensus_chr"]->m
 y[, i:=1:.N]

 y[, .(b=seq(start,end)), by=.(agp_chr, i)][, .N, key=.(agp_chr, b)][N > 1]->dups
 if(nrow(dups) > 0){
  stop(paste("Overlapping ranges specified: bin(s)",
        paste(dups[, paste(sep=":", agp_chr, b)], collapse=", "), "are contained in more than one inversion."))
 }

 y[, .(agp_chr=agp_chr[1], hic_bin = start:end, new_bin = end:start, new=T), by=i][, i := NULL][m, on=c("agp_chr", "hic_bin")][is.na(new), new := F]->m
 m[new == T, hic_bin := new_bin]
 m[new == T, consensus_orientation := consensus_orientation * -1]
 m[new == T & is.na(consensus_orientation), consensus_orientation := -1]
 m[, c("new", "new_bin", "agp_chr") := list(NULL, NULL, NULL)]
 make_agp(m, gap_size=hic_map$gap_size, species=species)->a

 list()->new
 new$hic_map <- m
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_multi_inversions <- copy(ranges)
 new
}

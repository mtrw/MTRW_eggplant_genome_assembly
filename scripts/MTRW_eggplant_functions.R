
library(igraph)
library(zoo)
library(stringi)
library(colorspace)



make_agp_file <- function( agp , file="plots/manual/manual_agp.csv") {
  write.csv(agp[scaffold!="gap"][order(agp_chr,agp_start),.(scaffold,orientation,idx=1:.N),by=.(agp_chr)],file,row.names = F,na = 0)
}


update_agp_from_file <- function(assembly,agp_file="plots/manual/agp_manual1.csv",gap_len=100,hic_nucfile="data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp_split.nuc.txt"){
  ass <- copy(assembly)
  agp_manual <- fread(agp_file) #cols: agp_chr and scaffold
  #just because it's stupid to get this from here when it's no doubt been messed up in editing. So delete and reload it directly from the time the SSs were made
  if("order" %in% colnames(agp_manual)){
    setorder(agp_manual,order)
    agp_manual[,order:=NULL]
  }
  if("index" %in% colnames(agp_manual)){
    agp_manual[,idx:=NULL]
  }
  if("flip" %in% colnames(agp_manual)){
    agp_manual[!is.na(flip) & !is.na(orientation),orientation:=-orientation]
  }
  agp_manual <- agp_manual[,.(agp_chr,scaffold,orientation)]
  agp_manual[is.na(orientation),orientation := 1]
  agp_manual[,idx:=1:.N,by=.(agp_chr)]
  fill_gaps <- copy(agp_manual)
  fill_gaps[,scaffold:="gap"]
  fill_gaps[,orientation:=1]
  fill_gaps <- fill_gaps[,.SD[1:(.N-1)],by=.(agp_chr)]
  agp_manual <- rbind(agp_manual,fill_gaps)[order(agp_chr,idx)]
  agp_manual <- assembly$info[,.(scaffold,scaffold_length = length,cM)][agp_manual,on="scaffold"]
  agp_manual[ scaffold == "gap" , scaffold_length := gap_len ]
  agp_manual[ , idx:=1:.N , by=.(agp_chr) ]
  
  #add AGP positions again (based on lengths)
  agp_manual[, agp_start := c(1,(cumsum(scaffold_length)[1:(.N-1)]+1) ) , by="agp_chr" ]
  agp_manual[, agp_end := cumsum(scaffold_length) , by="agp_chr" ]
  
  agp_manual[ scaffold != "gap" , hic_bin := 1:.N , by="agp_chr" ]
  agp_manual[, index := 1:nrow(agp_manual)  ]
  
  
  #
  #old_agp <- old_map$agp
  agp_manual[ , chr := gsub("[^0-9.-]","",agp_chr) %>% as.integer]
  agp_manual[,linkage_group := agp_chr]
  agp_manual[,gap := scaffold=="gap"]
  
  #old_agp_bed <- old_map$agp_bed
  agp_bed <- agp_manual[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), chr )]
  
  #old_chrlen <- old_map$chrlen
  chrlen <- agp_manual[scaffold != "gap" , .(length=max(agp_end))  , by=.(agp_chr,alphachr=chr,chr,truechr=!is.na(chr))]
  
  #old_hic_map <- old_map$hic_map
  hic_map <- NULL
  
  #old_hic_map_bin <- old_map$hic_map_bin
  hic_map_bin <- NULL
  
  
  binsize <- NULL
  gap_size <- gap_len #memory's sake
  max_cM_dist <- NULL
  
  manual_map <- list()
  manual_map$agp <- agp_manual
  manual_map$agp_bed <- agp_bed
  manual_map$chrlen <- chrlen
  manual_map$hic_map <- hic_map
  manual_map$binsize <- binsize
  manual_map$gap_size <- gap_size
  manual_map$max_cM_dist <- max_cM_dist
  ce("Manual map looks like this:")
  print(manual_map)
  ce("Adding AGP-wide analyses")
  manual_map <- add_psmol_fpairs_ep( assembly = assembly , hic_map = manual_map , nucfile = hic_nucfile ) 
  manual_map$hic_1Mb <- bin_hic_step(hic=manual_map$links, frags=manual_map$frags, binsize=2e5, chrlen=manual_map$chrlen, cores=21)
  manual_map$hic_1Mb$norm <- normalize_cis( manual_map$hic_1Mb, ncores=21, percentile=2 )
  manual_map$hic_1Mb$asymmetry <- find_inversions_ep( hic_map=manual_map,links=manual_map$hic_1Mb$norm , cores=30 )
  
  manual_map
}


read_fragdata_ep <- function( info , file ){
  fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
  fragbed[, length := end - start]
  fragbed[, start := start + 1]
  info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
  fragbed[, start := start - orig_start + 1]
  fragbed[, end := end - orig_start + 1]
  fragbed[, orig_start := NULL]
  fragbed[, orig_scaffold := NULL]
  fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
  list(bed=fragbed[], info=z[])
}
init_assembly_ep <- function( orig_agp , fai ){
  assembly <- list()
  
  info0 <- copy( orig_agp[,.(scaffold,orig_agp_chr)])[copy(fai),on=.(scaffold) ]
  info0[, orig_start := 1 ]
  info0[, orig_end := length ]
  info0[, orig_scaffold := scaffold ]
  
  assembly$info <- info0
  assembly[]
}



# anchor_scaffolds_ep <- function( assembly , nlinks_min=1 , r1rest_min=0.5 ){
#   a <- copy(assembly)
#   ce("Beginning length proportion anchored: ",(a$info[!is.na(orig_agp_chr),sum(length)]/a$info[,sum(length)]) %>% round(2) )
#   #attach chr info to link table
#   fp_id <- a$info[ !is.na(orig_agp_chr) , .(scaffold1=scaffold,scaff1_chr=orig_agp_chr) ][ a$fpairs[ scaffold1 != scaffold2 ] , on=.(scaffold1) ]
#   fp_id <- a$info[ !is.na(orig_agp_chr) , .(scaffold2=scaffold,scaff2_chr=orig_agp_chr) ][ fp_id , on=.(scaffold2) ]
#   
#   #counters
#   n_new_assigned <- 0
#   n_new_assigned_old <- -1
#   
#   #iterative ID loop  
#   #newly ID as many as possible: fp_id[only informative links][...]
#   while(n_new_assigned != n_new_assigned_old){
#     ce("N scaffolds assigned: ",(na <- a$info[!is.na(orig_agp_chr),.N]+n_new_assigned)," (",(na/nrow(a$info)) %>% round(digits = 2),"%)")
#     newids <- fp_id[is.na(scaff1_chr) & !is.na(scaff2_chr)][ , .(nlinks=.N) , by=.(scaffold1,scaff2_chr)][ order(nlinks) , {
#       if(length(nlinks)>1){
#         data.table(
#           n_scaff_links=sum(nlinks),
#           scaff_links_r12=nlinks[1]/(nlinks[1]+nlinks[2]),
#           scaff_links_r1rest=nlinks[1]/sum(nlinks),
#           best_chr=scaff2_chr[1]
#         )
#       } else {
#         data.table(
#           n_scaff_links=sum(nlinks),
#           scaff_links_r12=nlinks[1]/(nlinks[1]),
#           scaff_links_r1rest=nlinks[1]/sum(nlinks),
#           best_chr=scaff2_chr[1]
#         )
#       }
#     } , by=.(scaffold1) ]
#     #filter
#     newids_joiner <- newids[ nlinks_min>=nlinks_min & scaff_links_r1rest>=r1rest_min , .(scaffold1,scaffold2=scaffold1,best_chr) ]
#     
#     n_new_assigned_old <- n_new_assigned
#     n_new_assigned <- n_new_assigned+nrow(newids_joiner)
#     
#     #update fp_ids (scaff1 and then scaff2)
#     fp_id <- newids_joiner[,.(scaffold1,best_chr)][fp_id,on=.(scaffold1)]
#     fp_id[ !is.na(best_chr) , scaff1_chr:=best_chr ]
#     fp_id[ , best_chr := NULL ]
#     fp_id <- newids_joiner[,.(scaffold2,best_chr)][fp_id,on=.(scaffold2)]
#     fp_id[ !is.na(best_chr) , scaff2_chr:=best_chr ]
#     fp_id[ , best_chr := NULL ]
#     ce("\tPossible scaffolds to still assign: ",nu(fp_id[is.na(scaff1_chr)]$scaffold1))
#     #end loop
#   }
#   
#   ce("Convergence reached: updating chromosome assignments and returning updated assembly object.")
#   #update info
#   a$info <- fp_id[ !is.na(scaff1_chr) & !duplicated(scaffold1) , .(scaffold=scaffold1,best_chr=scaff1_chr)][a$info,on=.(scaffold) ]
#   a$info[,agp_chr:=orig_agp_chr]
#   a$info[ !is.na(best_chr) , agp_chr:=best_chr ]
#   a$info[ , best_chr := NULL ]
#   ce("Final length proportion anchored: ",(a$info[!is.na(orig_agp_chr),sum(length)]/a$info[,sum(length)]) %>% round(2) )
#   a
# }



anchor_scaffolds_gmap_ep <- function( assembly = Sm_v2 , nlinks_min=5 , r1rest_min=-Inf , r12_min=0.6 ){ #assumes GMAP is already done and used to anchor i.e. popseq_chr exists ... but Hi-C has not been used to expand this (hic_chr is absent)
  a <- copy(assembly)
  
  ce("Initial length proportion anchored: ",(a$info[!is.na(popseq_chr),sum(length)]/a$info[,sum(length)]) %>% round(2) )
  
  #attach chr info to link table
  fp_id <- a$info[ !is.na(popseq_chr) , .(scaffold1=scaffold,scaff1_chr=popseq_chr,scaff1_gmap_assigned=popseq_chr) ][ a$fpairs[ scaffold1 != scaffold2 ] , on=.(scaffold1) ]
  fp_id <- a$info[ !is.na(popseq_chr) , .(scaffold2=scaffold,scaff2_chr=popseq_chr,scaff2_gmap_assigned=popseq_chr) ][ fp_id , on=.(scaffold2) ]
  
  ce("Link counts: ")
  fp_id[,.N,by=.(scaff1_chr,scaff2_chr)] %>% dcast(scaff1_chr ~ scaff2_chr , value.var="N" )
  
  #counters
  n_new_assigned <- 0
  n_new_assigned_old <- -1
  
  #iterative ID loop
  #newly ID as many as possible: fp_id[only informative links][...]
  while(n_new_assigned != n_new_assigned_old){
    ce("N scaffolds assigned by Hi-C: ",(na <- a$info[!is.na(popseq_chr),.N]+n_new_assigned)," (",(na/nrow(a$info)) %>% round(digits = 2),"%)")
    newids <- fp_id[is.na(scaff1_chr) & !is.na(scaff2_chr)][ , .(nlinks=.N) , by=.(scaffold1,scaff2_chr)]
    setorder(newids,scaffold1,-nlinks)
    newids <- newids[ , {
      #if(scaffold1=="Super-Scaffold_46_a" ) {browser()}
      if(length(nlinks)>1){
        data.table(
          n_scaff_links=sum(nlinks),
          scaff_links_r12=nlinks[1]/(nlinks[1]+nlinks[2]),
          scaff_links_r1rest=nlinks[1]/sum(nlinks),
          best_hic_chr=scaff2_chr[1]
        )
      } else {
        data.table(
          n_scaff_links=sum(nlinks),
          scaff_links_r12=nlinks[1]/(nlinks[1]),
          scaff_links_r1rest=nlinks[1]/sum(nlinks),
          best_hic_chr=scaff2_chr[1]
        )
      }
    } , by=.(scaffold1) ]
    #filter
    newids_joiner <- newids[ n_scaff_links>=nlinks_min & scaff_links_r1rest>=r1rest_min & scaff_links_r12>=r12_min , .(scaffold1,scaffold2=scaffold1,best_hic_chr) ]
    
    n_new_assigned_old <- n_new_assigned
    n_new_assigned <- n_new_assigned+nrow(newids_joiner)
    
    
    #update fp_ids (scaff1 and then scaff2)
    fp_id <- newids_joiner[,.(scaffold1,best_hic_chr)][fp_id,on=.(scaffold1)]
    fp_id[ !is.na(best_hic_chr) , scaff1_chr:=best_hic_chr ]
    fp_id[ , best_hic_chr := NULL ]
    fp_id <- newids_joiner[,.(scaffold2,best_hic_chr)][fp_id,on=.(scaffold2)]
    fp_id[ !is.na(best_hic_chr) , scaff2_chr:=best_hic_chr ]
    fp_id[ , best_hic_chr := NULL ]
    ce("\tPossible scaffolds to still assign: ",nu(fp_id[is.na(scaff1_chr)]$scaffold1))
    #end loop
  }
  
  #now allow them to contradict their previous assignments to see if they converge
  n_reassigned <- 1
  
  while(n_reassigned > 0){
    ce("Reassigning ... ")
    newids <- fp_id[ !is.na(scaff1_chr) & !is.na(scaff2_chr)][ , .(nlinks=.N) , by=.(scaffold1,scaff1_chr,scaff2_chr)]
    setorder(newids,scaffold1,-nlinks)
    newids <- newids[ , {
      #if(scaffold1=="Super-Scaffold_46_a" ) {browser()}
      if(length(nlinks)>1){
        data.table(
          n_scaff_links=sum(nlinks),
          scaff_links_r12=nlinks[1]/(nlinks[1]+nlinks[2]),
          scaff_links_r1rest=nlinks[1]/sum(nlinks),
          old_chr=scaff1_chr[1],
          reassign_to=scaff2_chr[1]
        )
      } else {
        data.table(
          n_scaff_links=sum(nlinks),
          scaff_links_r12=nlinks[1]/(nlinks[1]),
          scaff_links_r1rest=nlinks[1]/sum(nlinks),
          old_chr=scaff1_chr[1],
          reassign_to=scaff2_chr[1]
        )
      }
    } , by=.(scaffold1) ]
    #filter
    newids_joiner <- newids[ (reassign_to != old_chr) & n_scaff_links>=nlinks_min & scaff_links_r1rest>=r1rest_min & scaff_links_r12>=r12_min , .(scaffold1,scaffold2=scaffold1,best_hic_chr=reassign_to,old_hic_chr=old_chr,scaff_links_r12,n_scaff_links) ]
    
    
    n_reassigned <- nrow(newids_joiner)
    ce("Reassigned ",n_reassigned," scaffolds this round.")
    
    #update fp_ids (scaff1 and then scaff2)
    
    fp_id <- newids_joiner[,.(scaffold1,best_hic_chr,scaff_links_r12_temp=scaff_links_r12,n_scaff_links_temp=n_scaff_links)][fp_id,on=.(scaffold1)]
    fp_id[ !is.na(best_hic_chr) , scaff1_chr:=best_hic_chr ]
    fp_id[ scaffold1 %in% newids_joiner$scaffold1 , scaff_links_r12 := scaff_links_r12_temp ]
    fp_id[ scaffold1 %in% newids_joiner$scaffold1 , n_scaff_links := n_scaff_links_temp ]
    fp_id[ , best_hic_chr := NULL ]
    fp_id[ , scaff_links_r12_temp := NULL ]
    fp_id[ , n_scaff_links_temp := NULL ]
    fp_id <- newids_joiner[,.(scaffold2,best_hic_chr)][fp_id,on=.(scaffold2)]
    fp_id[ !is.na(best_hic_chr) , scaff2_chr:=best_hic_chr ]
    fp_id[ , best_hic_chr := NULL ]
    #end loop
  }
  
  
  ce("Convergence reached. Making tough decisions, updating chromosome assignments, and returning assembly object")

  #update info
  new_info <- fp_id[ !is.na(scaff1_chr) & !duplicated(scaffold1) , .(scaffold=scaffold1,best_hic_chr=scaff1_chr,n_scaff_links,scaff_links_r12)][a$info,on=.(scaffold) ]
  
  #decision time
  new_info[ is.na(best_hic_chr) & !is.na(popseq_chr) , agp_chr := popseq_chr ] #forced move 1
  new_info[ !is.na(best_hic_chr) & is.na(popseq_chr) , agp_chr := best_hic_chr ] #forced move 2
  new_info[ best_hic_chr == popseq_chr , agp_chr := best_hic_chr ] #unambiguous
  new_info[ (best_hic_chr != popseq_chr) & (n_markers < 3) & (n_scaff_links < 10) ]$scaffold -> no_chr_scaffs #shit data all around
  new_info[ (best_hic_chr != popseq_chr) & (n_markers < 3) & (n_scaff_links < 10) , agp_chr := NA ] #shit data all around
  new_info[ (best_hic_chr != popseq_chr) & (n_markers >= 3) , agp_chr := popseq_chr ] #Gmap slightly more data than Hi-C
  new_info[ (best_hic_chr != popseq_chr) & (n_markers < 3) & (n_scaff_links >= 10) & (scaff_links_r12 > .8) ]$scaffold -> no_gmap_scaffs #Hi-C slightly more data than gmap
  new_info[ (best_hic_chr != popseq_chr) & (n_markers < 3) & (n_scaff_links >= 10) & (scaff_links_r12 > .8) , agp_chr := best_hic_chr ] #Hi-C slightly more data than gmap
  
  ce("Removing GMAP info for bad GMAP scaffs: ")
  print(no_chr_scaffs)
  print(no_gmap_scaffs)
  
  #remove gmap info for scaffs where that info is deemed unreliable or overwritten by Hi-C (as above)
  new_info[scaffold %in% no_chr_scaffs , cM := NA ]
  new_info[scaffold %in% no_chr_scaffs , linkage_group := NA ]
  new_info[scaffold %in% no_chr_scaffs , popseq_chr := NA ]
  new_info[scaffold %in% no_gmap_scaffs , cM := NA ]
  new_info[scaffold %in% no_gmap_scaffs , linkage_group := NA ]
  new_info[scaffold %in% no_gmap_scaffs , popseq_chr := NA ]
  #browser()
  a$fragdata$info[ , cM := NULL ]
  a$fragdata$info[ , linkage_group := NULL ]
  a$fragdata$info[ , popseq_chr := NULL ]
  new_fragdata_info <- new_info[, .(agp_chr,scaffold,cM,linkage_group,popseq_chr) ][a$fragdata$info,on="scaffold"]
  
  a$info <- new_info
  a$fragdata$info <- new_fragdata_info
  
  ce("Final length proportion anchored: ",(a$info[!is.na(agp_chr),sum(length)]/a$info[,sum(length)]) %>% round(2) )
  a
}

chrNames <- function(agp=NULL,species=NULL) {
  data.table( alphachr=paste0("linkage_group_",1:15) , chr=1:15L )
}


find_inversions_ep <- function(hic_map, links, species="wheat", cores=1, winsize=15, maxdist=1e8, threshold=40, factor=100){
  cat("In find inversions routine \n")
  
  unique(links$chr1) -> chrs
  
  links[chr1 %in% chrs]->zz
  zz[, l := log(factor*nlinks_norm/sum(nlinks_norm)*1e6)]->zz
  hic_map$chrlen[, .(chr1=chr, length)][zz, on="chr1"]->zz
  zz[dist <= bin1 & dist <= bin2 & (length - bin1) >= dist & (length - bin2) >= dist]->zz
  zz[l >= 0 & dist <= maxdist, .(r = sum(sign(bin1 - bin2) * l)), key=.(chr=chr1, bin=bin1)]->w
  chrNames_ep(agp=T)[w, on="chr"]->w
  
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


add_psmol_fpairs_ep <- function(assembly, hic_map, nucfile ){
  
  assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
  
  hic_map$agp[agp_chr != "chrUn", .(chr=chr, scaffold, orientation=orientation, agp_start, agp_end)]->a
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




















chrNames_ep <- function(agp=F, species="eggplant") {
  if(species == "wheat"){
    data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
  } else if (species == "barley"){
    data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
  }
  else if (species == "rye"){
    data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
  }
  else if (species == "eggplant"){
    data.table(alphachr=apply(expand.grid("E",1:12, stringsAsFactors=F), 1, function(x) paste0(gsub("\\s*","",x), collapse="") ), chr=1:12)->z
  }
  if(agp){
    rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
  }
  z[]
}
#


chrNames_gmap_ep <- function(agp=F) {
  data.table(linkage_group=c("E1","E2A","E2B","E3","E4","E5","E6","E7","E8A","E8B","E9","E10A","E10B","E11","E12"), chr=1:15L )
}








bin_hic_step<-function(hic, frags, binsize, step=NULL, chrlen, cores=1){
  
  chrs <- unique(chrlen$chr)
  
  if(is.null(step)){
    step <- binsize
  }
  
  options(scipen=20)
  chrlen[chr %in% chrs, .(bin=seq(0, length, step)), key=chr][, end := bin + binsize - 1][, start := bin] -> bins
  bins[, id := paste0(chr, ":", bin)]
  
  hic[chr1 %in% chrs & chr2 %in% chrs] -> z
  
  if(step == binsize){
    z[, bin1 := as.integer(start1 %/% binsize * binsize)]
    z[, bin2 := as.integer(start2 %/% binsize * binsize)]
    z[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]->z
    z[, id1 := paste(sep=":", chr1, bin1)]
    z[, id2 := paste(sep=":", chr2, bin2)]
    z->mat
    
    frags[chr %in% chrs]->f
    f[, bin := as.integer(start %/% binsize * binsize)]
    f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
          gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
    f[, id := paste(sep=":", chr, bin)]
  } else {
    copy(bins)->b
    b[, idx := 1:.N]
    b[, .(win=seq(bin, bin+binsize-1, step)),.(chr, bin)]->b
    
    z[chr1 == chr2]->z
    z[, win1:=as.integer(start1 %/% step * step)]
    z[, win2:=as.integer(start2 %/% step * step)]
    
    rbindlist(mclapply(mc.cores=cores, mc.preschedule=F, chrs, function(i){
      z[chr1 == i]->y
      y[, .(nlinks=.N), keyby=.(chr1,win1,chr2,win2)]->y
      setnames(copy(b), paste0(names(b), 1))[y, on=c("chr1", "win1"), allow.cartesian=T]->y
      setnames(copy(b), paste0(names(b), 2))[y, on=c("chr2", "win2"), allow.cartesian=T]->y
      y[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]
    }))->y
    y[, id1:=stri_c(sep=":", chr1, bin1)]
    y[, id2:=stri_c(sep=":", chr2, bin2)]
    y[, dist := abs(bin1 - bin2)]
    y->mat
    
    copy(frags)->f
    f[, win := as.integer(start %/% step * step)]
    b[f, on=c("chr", "win"), nomatch=0, allow.cartesian=T]->f
    f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
          gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
    f[, id :=paste(sep=":", chr, bin)]
  }
  
  list(bins=f[], mat=mat[])
}












normalize_cis<-function(binhic , ncores=1 , chrs=NULL , percentile=0 , model_inc_distance=FALSE , model_inc_mapability=FALSE , KRnormalise=FALSE , distance_cancel=F){
  if( is.null(chrs) ){
    chrs <- unique(binhic$bins$chr)
  }
  
  mf <- binhic$bins[chr %in% chrs]
  ab <- binhic$mat
  
  if(percentile > 0){
    mf[eff_length >= quantile(eff_length, 0:100/100)[percentile + 1]]->mf
  }
  
  #dev: i <- 1
  rbindlist(mclapply(mc.cores=ncores, chrs, function(i){
    ab[chr1 == chr2 & chr1 == i & id1 %in% mf$id & id2 %in% mf$id]->abf
    dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat
    
    u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
    setkey(mf, id)[colnames(u)]->v
    
    #add distance matrix to add to model? seems pretty fucking odd to pool them all together!
    if (model_inc_distance) {
      abf[ , dist := abs(bin1-bin2) ]
      allbins <- unique(c(abf$bin1,abf$bin2))
      diffgrid <- grid_ply(allbins,allbins,difff)
      dimnames(diffgrid) <- lapply(dimnames(diffgrid) , function(x) paste0(i,":",x) )
      diffgrid <- diffgrid[sort(rownames(diffgrid)),sort(colnames(diffgrid))]
      dim(diffgrid)
      if (sum(colnames(diffgrid)!=colnames(u)) + sum(rownames(diffgrid)!=rownames(u)) == 0){
        u_d <- diffgrid
      } else {
        stop("diff_matrix peculiar, i.e. wrong\n")
      }
    }
    
    if (model_inc_distance) {
      normalize_mat_cis(u,v,u_d=u_d, plots=FALSE,model_inc_mapability=model_inc_mapability,KRnormalise=KRnormalise,distance_cancel=distance_cancel)
    } else {
      normalize_mat_cis(u,v,u_d=NULL,plots=FALSE,model_inc_mapability=model_inc_mapability,KRnormalise=KRnormalise,distance_cancel=distance_cancel)
    }
  }))->nhic
  mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
  mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
  nhic[, dist := abs(bin2 - bin1)][]
}



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################



normalize_mat_cis<-function( u , v , u_d=NULL , plots=FALSE , KRnormalise=FALSE , model_inc_mapability=F , distance_cancel=F ){
  u_vec<-u[upper.tri(u,diag=F)]
  
  #get cov matrix
  len_m<-as.matrix(log(v[, eff_length]%o%v[, eff_length]))
  gcc_m<-as.matrix(log(v[, gc]%o%v[, gc]))
  if( model_inc_mapability & !is.null(v$nonmappability_scaled ) ){
    map_m<-as.matrix(log(v[, nonmappability_scaled]%o%v[, nonmappability_scaled ]))
  } else {
    map_m<-as.matrix(log(v[, cov]%o%v[, cov]))
  }
  
  #centralize cov matrix of enz, gcc
  len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
  gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))
  if(  model_inc_mapability & !is.null(v$nonmappability_scaled) ){ #if it has actual values in it. If not, it's all zeroes and normalisation gives NaN
    map_m<-(map_m-mean(c(map_m)))/sd(c(map_m))
  }
  
  
  #change matrix into vector
  len_vec<-len_m[upper.tri(len_m,diag=F)]# %>% z_transform
  gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]# %>% z_transform
  map_vec<-map_m[upper.tri(map_m,diag=F)]# %>% z_transform
  if ( !is.null(u_d) ){ #if we're using a distance matrix
    logdist_mat <- log(u_d+1)
    logdist_mat<-(logdist_mat-mean(c(logdist_mat)))/sd(c(logdist_mat))
    logdist_vec <- logdist_mat[upper.tri(logdist_mat,diag=FALSE)]
    #d_vec<-u_d[upper.tri(u_d,diag=F)] #kill me
  } else {
    logdist_vec <- rep(0,length(len_vec))
  }
  #logdist_vec<-d_vec %>% log %>% `-` %>% z_transform #hhhhere these transforms should be done at the MATRIX stage
  
  #remember if it has been calculated, d_vec also exists in this set
  
  d <- data.table(eff_length=len_vec,links=u_vec,gcc=gcc_vec,mappability=map_vec,log_dist=logdist_vec)
  if ( !is.null(v$nonmappability_scaled) ){
    t <- as.matrix(log(v[, nonmappability_scaled]%o%v[, nonmappability_scaled ]))
    t<-(t-mean(c(t)))/sd(c(t))
    d[ , mappability_plot := t[upper.tri(t,diag=F)] ]
  }
  
  
  #fit Poisson regression: u~len+gcc+offset(map)
  if( model_inc_mapability & !is.null(v$nonmappability_scaled) & !is.null(u_d) ){
    fit<-glm(data=d,links~eff_length+gcc+mappability+log_dist,family="poisson")
    d[ , prediction := predict(fit , d[, .(eff_length,gcc,mappability,log_dist) ])]
    #for presentation
    if(plots) ggplot(d, aes(x=eff_length,y=exp(prediction)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,colour="red",alpha=1)
    #for pres
    if(plots) ggplot(d[log_dist < -2.5], aes(x=eff_length,y=exp(prediction),colour=as.factor(log_dist)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,alpha=1)
    
    if(plots) ggplot(d[log_dist < -2.5], aes(x=mappability,y=exp(prediction),colour=as.factor(log_dist)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,alpha=1)
    if(plots) ggplot(d, aes( x=links,y=links-exp(prediction) ) ) + geom_point(size=.5)
    #for presentation
    if(plots) ggplot(d, aes( x=links,y=exp(prediction) ) ) + geom_point(size=.5)
    coeff <- fit$coefficients
    res <- round(u/exp(coeff["eff_length"]*len_m + coeff["gcc"]*gcc_m + coeff["mappability"]*map_m), 4)
    if( distance_cancel ) {res <- round(u/exp(coeff["eff_length"]*len_m + coeff["gcc"]*gcc_m + coeff["mappability"]*map_m + coeff["log_dist"]*logdist_mat), 4) }
    #res3 <- round(u/exp(coeff["eff_length"]*len_m + coeff["gcc"]*gcc_m + coeff["mappability"]*map_m + coeff["log_dist"]*log(u_d)), 4)
  } else if ( model_inc_mapability & !is.null(v$nonmappability_scaled) ) { #model and adjust for mappability
    fit<-glm(data=d,links~eff_length+gcc+mappability,family="poisson")
    d[ , prediction := predict(fit , d[, .(eff_length,gcc,mappability) ])]
    if(plots) ggplot(d, aes(x=mappability,y=exp(prediction)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,colour="red",alpha=.1)
    if(plots) ggplot(d, aes( x=links,y=links-exp(prediction) ) ) + geom_point(size=.5)
    coeff <- fit$coefficients
  } else { #what we originally were doing (HiCNorm method without mappability)
    fit<-glm( data=d , links~eff_length+gcc+offset(mappability) , family="poisson" )
    coeff <- fit$coefficients
    d[ , prediction := coeff["(Intercept)"] + coeff["eff_length"]*eff_length + coeff["gcc"]*gcc + mappability ]
    #for pres
    if(plots) ggplot(d, aes(x=eff_length,y=exp(prediction)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,colour="red",alpha=.1)
    #for pres
    if(plots) ggplot(d, aes( x=links,y=exp(prediction) ) ) + geom_point(size=.5)
    #for pres
    if(plots) ggplot(d, aes(x=mappability_plot,y=exp(prediction)) ) + geom_point(size=.5) + geom_point(aes(y=links),size=.1,colour="red",alpha=.1)
    res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
  }
  
  if(KRnormalise == TRUE) { res <- KRnorm(res)}
  
  data.table(id1=colnames(res), res)->res
  melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
  res[nlinks_norm > 0]
}
#



plot_asymmetry_agp <- function(agp,asymmetry){
  ggplot(asymmetry$ratio,aes(x=bin+5e5,y=r)) + 
    geom_point(size=.1) +
    geom_vline(data=agp,aes(xintercept=agp_start),size=.1,alpha=.2) + 
    geom_text(data=agp , aes(x=agp_start,y=-100,label=scaffold),angle=90,size=1,colour="black" ) +
    facet_grid(chr~.)
}












plot_contact_matrices_ep <-function(hicmap){
  ggplot(hicmap$hic_1Mb$norm[chr1==chr2],aes(x=bin1,y=bin2)) +
    geom_tile(aes(fill = log10(nlinks_norm))) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme(panel.background = element_blank()) + 
    theme(legend.position="none") +
    facet_grid( chr1 ~ . ) + 
    geom_vline(data=hicmap$agp[,.(agp_start,chr1=chr)],aes(xintercept=agp_start),size=.05,alpha=.1) +
    geom_text(data=hicmap$agp[,.(agp_start,chr1=chr,scaffold_length,scaffold)] , aes(x=agp_start,y=1,label=scaffold),angle=90,size=.1,colour="black")
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





plot_cmp_agps <- function(agp_list,ask="Q",plotpos=F,colour_table=NULL){
  al <- copy(agp_list)
  data <- lapply(1:length(al) , function(x){
    z <- al[[x]][,.(scaffold,agp_start,agp_end,chr,orientation)]
    z[, agp:=x ]
    z
  }) %>% rbindlist()
  data <- data[!is.na(chr) & scaffold != "gap" ]
  data[ , new:=FALSE ]
  data[ , delete:=FALSE ]
  data[is.na(orientation) , orientation:=1]
  #find which drop out in next
  l_ply(1:(length(al)-1) , function(x) {
    data[agp==x & ( !scaffold %in% data[agp==(x+1)]$scaffold ) , delete := TRUE]
  })
  #which come in new
  l_ply(2:length(al) , function(x) {
    data[ agp==x & ( !scaffold %in% data[agp==(x-1)]$scaffold ) , new := TRUE ]
    b <- data[ agp==(x-1) , .( scaffold , p1=pmean(agp_start,agp_end) )  , by=.(chr) ]
    n <- data[ agp==x ,     .( scaffold , p2=pmean(agp_start,agp_end) )  , by=.(chr) ]
    cmp <- b[n,on=c("chr","scaffold")]
    #something odd here ... check whether it's in cmp or before. lines of points with ...
    cmp[ (p1+p2)<40e6 ][order(p1)]
    
    
    if(plotpos){
      par(mfrow = c(ceiling(nu(data$chr)/4),4),oma=rep(0,4),mar=rep(2.5,4))
      cmp[ , plot( p2 ~ p1 , .SD , pch=20 , xlab="p0" ) , by=.(chr)]
      { invisible(ask<-readline(prompt="Press [enter] to kill plot and continue. ")) }
      dev.off()
    }
    
  })
  
  nagp <- nu(data$agp)
  
  data[ , vpos := ((chr-1)*nagp) + ( ((nagp*.8)/(nagp-1)) * (agp-1) ) ]
  vinc_peragp <- ((nagp*.8)/(nagp-1)) / max(data$vpos) *  max(data$agp_end)
  data[ , vpos := (vpos/max(vpos)) * max(agp_end) ]
  data[ , vpos := max(vpos)-vpos]
  
  #add colours with colorspace
  if(is.null(colour_table)){
    if(nu(data$agp)>12) {stop("More than 12 AGPs supplied, please provide your own colour table ... ")}
    schemes <- (expand.grid( c("Greens","Blues","Reds","Purples") , c("",2,3) ) %>% apply(1,function(x) paste(x,collapse = " ")) %>% sub(" $","",.))
    colour_table <- data[ , .(colour=sequential_hcl(3,palette=schemes[agp],alpha=.8)[c(1,2)]) , by=.(agp) ]
  }
  
  data[,colour:=character()]
  setkey(data,agp,chr,agp_start)
  l_ply(u(data$agp),function(a){
    cols <- colour_table[agp==a]$colour
    l_ply(u(data$chr) , function(c){
      data[ agp==a & chr==c , colour := rep(cols,length.out=.N) ][]
    })
  })
  
  
  hrange <- range(data$agp_start,data$agp_end)
  vrange <- range(data$vpos)
  
  plot(
    NULL,
    xlim=hrange,
    ylim=vrange,
    ylab=NA,
    xlab="AGP position",
    yaxt="n"
  )
  
  l_ply(1:(nagp-1),function(x){
    data[agp==x | agp==(x+1) , {
      block <- .SD[order(agp)]
      #print(block)
      if( .N==1 & all(delete) ){
        #ce("Deleter at ",scaffold)
        polygon(x=c(agp_start,agp_end,mean(agp_start,agp_end)),y=c(vpos,vpos,vpos-(vinc_peragp/4)),border="red")
      } else if ( .N==1 & all(!delete) ) {
        ce("Adder at ",scaffold)
        polygon(x=c(agp_start,agp_end,mean(agp_start,agp_end)),y=c(vpos,vpos,vpos+(vinc_peragp/4)),border="darkgreen")
      } else if ( .N != 2 | nu(agp) != 2 ) {
        stop("Something's screwed up at",scaffold)
      } else {
        #ce("We have a block at ",scaffold)
        if(ask != "Q"){ invisible(ask<-readline(prompt="Press [enter] to continue, \"Q\" to stop asking ")) }
        if(sum(orientation)==0){
          #bowtie
          polygon(c(agp_start[1],agp_end[1],agp_start[2],agp_end[2]),c(vpos[1],vpos[1],vpos[2],vpos[2]),col=colour[1],lwd=.5,border="#1B1B1BCC")
        } else {
          #rectangle
          polygon(c(agp_start[1],agp_end[1],agp_end[2],agp_start[2]),c(vpos[1],vpos[1],vpos[2],vpos[2]),col=colour[1],lwd=.5,border="#1B1B1BCC")
        }
      }
    } , by=.(chr,scaffold) ]
  })
  data[ , abline(h=vpos[1]) , by=.(vpos) ]
  middle <- mean(data$agp_end)
  data[  , .SD[vpos==min(vpos) , text(x=middle , y=vpos[1]-(vinc_peragp/8), labels=paste("Chromosome",chr[1]) , cex=.5 ) ] , by=.(chr) ]
  #
}
#








add_hic_cov<-function(assembly, scaffolds=NULL, binsize=1e3, binsize2=1e6, minNbin=50, innerDist=1e5, cores=1){
  
  info<-assembly$info
  
  if("mr" %in% colnames(info) | "mri" %in% colnames(info)){
    stop("assembly$info already has mr and/or mri columns; aborting.")
  }
  
  fpairs<-assembly$fpairs
  
  if(is.null(scaffolds)){
    scaffolds <- info$scaffold
    null=T
  } else {
    info[scaffold %in% scaffolds]->info
    fpairs[scaffold1 %in% scaffolds]->fpairs
    null=F
  }
  
  fpairs[scaffold1 == scaffold2 & pos1 < pos2][, .(scaffold = scaffold1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
  f[bin2 - bin1 > 2*binsize]->f
  f[, i := 1:.N]
  f[, b := paste0(scaffold, ":", bin1 %/% binsize2)]
  setkey(f, b)
  
  rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
    f[j][, .(scaffold=scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
  }))->ff
  
  if(nrow(ff) > 0){
    ff[, .(n=sum(n)), key=.(scaffold, bin)]->ff
    info[, .(scaffold, length)][ff, on="scaffold"]->ff
    ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
    ff[, nbin := .N, key="scaffold"]
    ff[, mn := mean(n), key=d]
    ff[, r := log2(n/mn)]
    ff[nbin >= minNbin, .(mr=suppressWarnings(min(r))), key=scaffold][order(mr)]->z
    ff[nbin > minNbin & d >= innerDist, .(mri=suppressWarnings(min(r))), key=scaffold]->zi
    z[ff, on="scaffold"]->ff
    zi[ff, on="scaffold"]->ff
    z[info, on="scaffold"]->info_mr
    zi[info_mr, on="scaffold"]->info_mr
  } else {
    copy(info) -> info_mr
    info_mr[, c("mri", "mr") := list(NA, NA)]
  }
  
  if(null){
    assembly$info=info_mr
    assembly$cov=ff
    
    assembly$binsize <- binsize
    assembly$minNbin <- minNbin
    assembly$innerDist <- innerDist
    
    assembly
  } else {
    list(info=info_mr, cov=ff)
  }
}





















# Create diagnostics plot for putative chimeras
plot_chimeras<-function(scaffolds, assembly, breaks=NULL, file, mindist=0, cores=20, species="wheat"){
  
  plot_popseq_carma_tcc<-function(scaffold, breaks=NULL, page, info, cssaln, tcc_pos, span, molcov, species){
    chrNames(species=species)->wheatchr
    
    i <- scaffold
    
    if(is.null(tcc_pos) || nrow(tcc_pos) == 0){
      hic <- F
    } else {
      hic <- T
    }
    
    if(is.null(molcov) || nrow(molcov) == 0){
      tenex <- F
    } else {
      tenex <- T
    }
    
    nrow <- 1
    if(hic){
      nrow <- nrow + 1
    }
    if(tenex){
      nrow <- nrow + 1
    }
    
    par(mfrow=c(nrow,3))
    
    par(oma=c(0,0,3,0))
    if(hic){
      par(mar=c(1,4,4,1))
    } else {
      par(mar=c(4,4,4,1))
    }
    par(cex=1)
    
    if(!is.null(breaks)){
      br <- breaks[i, on="scaffold", br]/1e6
    } else if (hic){
      span[i, on="scaffold"][order(r)][1, bin]/1e6->br
    } else {
      br <- NULL
    }
    
    if(hic){
      info[scaffold == i, paste0(", bad Hi-C: ", bad_hic)] -> badhic
      xlab<-""
    } else {
      badhic <- ""
      xlab<-"position in scaffold (Mb)"
    }
    
    cssaln[i, on="scaffold"]->z
    l=info[i, on="scaffold"][, length]/1e6
    
    if(species %in% c("rye", "lolium", "barley")){
      ymin <- 8
    } else {
      ymin <- 24
    }
    
    plot(xlim=c(0,l), 0, ylim=c(ymin,1), bty='l', col=0, yaxt='n', pch=20, xlab=xlab, ylab="POPSEQ chromosome")
    if(!is.null(br)){
      abline(v=br, lwd=3, col="blue")
    }
    title("POPSEQ", line=1)
    z[, points(pos/1e6, popseq_chr, pch=20)]
    axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
    
    info[scaffold == i, title(outer=T, line=0, paste("page: ", page, ", ", scaffold, sep="", " (", round(l,1), " Mb),\n",  
                                                     "bad POPSEQ: ", bad_popseq, ", bad CARMA: ", bad_sorted, badhic))]
    
    w<-info[i, on="scaffold"]$popseq_chr
    if(is.na(w)){
      plot(xlab='', ylab='', type='n', 0, axes=F)
    } else {
      z[popseq_chr == w][, plot(pch=20, xlab=xlab, ylab="genetic position (cM)", xlim=c(0,l), bty='l', las=1, pos/1e6, popseq_cM)]
      title(main="Genetic positions of CSS contigs\nfrom major chromosome")
      if(hic){
        abline(v=br, lwd=3, col="blue")
      }
    } 
    
    if(species != 'lolium'){
      plot(xlim=c(0,l), 0, bty='l', ylim=c(ymin,1), yaxt='n', pch=20, col=0, xlab=xlab, ylab="flow-sorting chromosome")
      title("flow-sorting CARMA", line=1)
      if(!is.null(br)){
        abline(v=br, lwd=3, col="blue")
      }
      if(species == "wheat"){
        z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "3B" | sorted_arm == "S", 1, 2), pch=20)]
        legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c("short arm / 3B", "long arm"))
      }
      if(species == "rye"){
        z[, points(pos/1e6, sorted_chr, col=1, pch=20)]
      }
      if(species == "barley"){
        z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "1H" | sorted_arm == "S", 1, 2), pch=20)]
        legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c("short arm / 1H", "long arm"))
      }
      axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
    } else {
      plot(0, type='n', xlab="", ylab="", axes=F)
    }
    
    if(hic){
      par(mar=c(4,4,4,1))
      plot(0, xlim=c(0,l), col=0, bty='l', ylab='Hi-C chromosome', yaxt='n', ylim=c(ymin,1), xlab="position in scaffold (Mb)")
      title("Interchromosomal Hi-C links", line=1)
      tcc_pos[i, on="scaffold1", points(pos1/1e6, col="#00000003", chr2, pch=20)]
      abline(v=br, lwd=3, col="blue")
      axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
      
      span[i, on="scaffold"]->w
      if(nrow(w[!is.na(n)]) == 0){
        plot(xlab='', ylab='', type='n', 0, axes=F)
      } else {
        w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
        title("Intrascaffold physical Hi-C coverage", line=1)
        abline(v=br, lwd=3, col="blue")
        w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
      }
      
      if(nrow(w[!is.na(n)]) == 0){
        plot(xlab='', ylab='', type='n', 0, axes=F)
      } else {
        w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
        title("Hi-C expected vs. observed coverage", line=1)
        abline(v=br, lwd=3, col="blue")
        w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
      }
    }
    
    if(tenex){
      molcov[i, on="scaffold"]->w
      if(nrow(w) == 0){
        plot(xlab='', ylab='', type='n', 0, axes=F)
      } else {
        w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
        title("10X molecule coverage", line=1)
        if(!is.null(br)){
          abline(v=br, lwd=3, col="blue")
        }
        w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
      }
      
      if(nrow(w) == 0){
        plot(xlab='', ylab='', type='n', 0, axes=F)
      } else {
        w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
        title("10X expected vs. observed coverage", line=1)
        if(!is.null(br)){
          abline(v=br, lwd=3, col="blue")
        }
        w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
      }
    }
  }
  
  info<-assembly$info
  ff<-assembly$cov
  tcc_pos <- assembly$fpairs
  cssaln <- assembly$cssaln
  molcov <- assembly$molecule_cov
  
  width <- 2500
  height <- 1000
  res <- 150
  
  if(!is.null(tcc_pos) && nrow(tcc_pos) > 0){
    height <- height + 1000
  }
  if(!is.null(molcov) && nrow(molcov) > 0){
    height <- height + 1000
  }
  
  scaffolds <- scaffolds[!duplicated(scaffold), .(scaffold)]
  
  bad <- copy(scaffolds)
  out <- file
  bad[, f:=tempfile(fileext=".png"), by=scaffold]
  bad[, p:=tempfile(fileext=".pdf"), by=scaffold]
  mclapply(mc.cores=cores, 1:nrow(bad), function(i){
    file <- bad[i, f]
    pdf <- bad[i, p]
    s<-bad[i,scaffold]
    cat(paste0(i, " ", bad[i, scaffold], "\n"))
    png(file, height=height, res=res, width=width)
    plot_popseq_carma_tcc(s, breaks=breaks, page=i, info, cssaln, tcc_pos, ff[d >= mindist], molcov, species)
    dev.off()
    system(paste("convert", file, pdf))
    unlink(file)
  })
  system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", out, " ", paste(bad$p, collapse=' ')))
  unlink(bad$p)
}

plot_genetic_map_vs_agp <- function(assembly,hicagp,transpose=FALSE,map=NULL){
  m <- map_to_agp_coords(assembly,hicagp,transpose,map)
  ggplot( data=m , aes(x = marker_agp_pos , y = marker_popseq_cM , colour = as.factor(marker_popseq_lg_chr))) + 
    geom_point(data=m , size=.2)  + 
    geom_vline(data=hicagp , mapping=aes(xintercept = agp_start),alpha = .2 , colour="black") + 
    facet_grid(agp_chr~.) +
    geom_text(data=hicagp , aes(x=agp_start,y=100,label=ifelse(scaffold_length > 500 , sub(".*_","",scaffold) , "")),angle=90,size=3,colour="black")
}


plot_genetic_map_vs_agp_ep <- function(markers,agp,return_data=F){
  nchr <- length(unique(agp$chr))
  mpos <- agp[,.(scaffold,agp_chr=chr,agp_start,agp_end,orientation)][markers,on="scaffold"]
  mpos[ , colour := replace_levels_with_colours(linkage_group,fun="qualitative_hcl",palette = "Dark 3") ]
  mpos[,pos_in_agp := ifelse( orientation > 0 , (agp_start + pos_in_scaff) - 1 , (agp_end - pos_in_scaff) + 1 )]
  if(return_data){return(mpos)}
  dev.off()
  .old_par <- par()
  par(mfrow=(c(4,1)),mar=c(2.7,2,2,2))
  for(ch in sort(unique(agp$chr))){
    null_plot(
      x = c(1,max(agp$agp_end)),
      y = range(markers$cM,na.rm = T),
      xlab=NA,
      ylab=NA
    )
    title(
      xlab = paste0("Position chromosome ",ch),
      ylab = "cM",
      line = 2
    )
    abline( v = agp[scaffold!="gap"][chr==ch,c(agp_start,agp_end)] , lwd=.2 , col="#00000055" )
    
    points(
      x = mpos[agp_chr==ch]$pos_in_agp,
      y = mpos[agp_chr==ch]$cM,
      col = mpos[agp_chr==ch]$colour,
      cex = .5,
      pch = 20
    )
    if(! ch %% 4 & ch != nchr){
      wait(message="For next four chromosomes, please press [enter].")
    }
  }
  par(.old_par)
}

plot_many_maps_vs_agp_ep <- function(markers,agp,return_data=F){
  .old_par <- par()
  nchr <- length(unique(agp$chr))
  nmaps <- length(unique(markers$map))
  nlgs <- length(unique(markers$linkage_group))
  mpos <- agp[,.(scaffold,agp_chr=chr,agp_start,agp_end,orientation)][markers,on="scaffold"]
  maplist <- sort(unique(mpos$map))
  mpos[ , colour := replace_levels_with_colours(linkage_group,fun="qualitative_hcl",palette = "Dark 3") ]
  colour_table <- mpos[,.N,by=.(colour,linkage_group)]
  setkey(colour_table,"linkage_group")
  mpos[,pos_in_agp := ifelse( orientation > 0 , (agp_start + pos_in_scaff) - 1 , (agp_end - pos_in_scaff) + 1 )]
  if(return_data){return(mpos)}
  dev.off()
  par(mfrow=(c(5,nmaps)),mar=c(2.7,2,2,2))
  layout(matrix(c(rep(1,nmaps),2:((nmaps*4)+1)),ncol=nmaps,byrow=T))
  printcols <- TRUE
  for(ch in sort(unique(agp$chr))){
    if(printcols){
      null_plot(1:nlgs,1,xaxt="n",yaxt="n")
      points(1:nlgs,rep(1,nlgs),col=colour_table$colour,cex=2,pch=20)
      text(1:nlgs,rep(0.8,nmaps),labels=colour_table$linkage_group,cex=.6)
      printcols <- 0
    }
    for(mp in maplist){
      rowdat <- mpos[map==mp]
      null_plot(
        x = c(1,max(agp$agp_end)),
        y = range(markers$cM,na.rm = T),
        xlab=NA,
        ylab=NA
      )
      title(
        xlab = paste0("Position chromosome ",ch," (map ",mp,")"),
        ylab = "cM",
        line = 2
      )
      abline( v = agp[scaffold!="gap"][chr==ch,c(agp_start,agp_end)] , lwd=.2 , col="#00000015" )
      par(srt=-90)
      text(
        x = agp[scaffold!="gap"][chr==ch,agp_start] ,
        y = agp[scaffold!="gap"][chr==ch,rep(mean(cM,na.rm=TRUE),.N)] ,
        labels = agp[scaffold!="gap"][chr==ch,scaffold] ,
        cex = 0.02,
        adj = c(0,1)
      )
      par(srt=0)
      points(
        x = rowdat[agp_chr==ch]$pos_in_agp,
        y = rowdat[agp_chr==ch]$cM,
        col = rowdat[agp_chr==ch]$colour,
        cex = .5,
        pch = 20
      )
      if(mp==last(maplist) & ! ch %% 4 & ch != nchr){
        wait(message="For next four chromosomes, please press [enter].")
        printcols <- 1
      }
    }
  }
  par(.old_par)
}

#Load packages
setwd("/filer-dg/agruppen/seq_shared/MTRW_eggplant_assembly/")
source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R")
source("scripts/MTRW_eggplant_functions.R")
source("scripts/hic_mapping_functions_ep.R")


#Get scaffold fai
fai <- fread("data/scaffolds/Scaffold_to_HIC.fa.fai",col.names = c("scaffold","length") , select = 1:2 )

#Get original AGP
orig_agp <- fread("data/Scaffolds_superscaffolds_eggplant_to_SMELV3.agp",col.names=c("scaffold","agp_orientation","orig_agp_chr","agp_start"))
orig_agp[ , agp_start:=agp_start+1 ] #make 1-based
orig_agp[ , orig_agp_chr:=sub("SMEL3Ch","",orig_agp_chr) %>% as.integer ]

#check AGP is all real chrs
orig_agp$scaffold %prop_in% fai$scaffold

#initialise assembly object
Sm_v1 <- init_assembly_ep( orig_agp=orig_agp , fai=fai )

#add hic data objects
Sm_v1$fpairs <- fread("zcat mapping/MTRW_eggplant_hic_fragment_pairs.tsv.gz", header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))
Sm_v1$fragdata <- read_fragdata_ep( info=Sm_v1$info , file="data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp.bed" )


#hic QC
dist_data <- ldply(c("mapping/*_length_dist_PE.txt","mapping/*_length_dist.txt"),function(d){
  a <- Sys.glob(d)
  dt <- fread(a,col.names = c("count","length"))
  dt[ , PE:=grepl("_PE",d) ]
  dt
}) %>% as.data.table

# ggplot( dist_data , aes(x=length,y=count,colour=PE)) +
#   geom_line() -> p; p
# ggsave(filename = file.path("plots","hic_length_dist.pdf") , plot = p , device = "pdf" , units = "in" , width = 6 , height = 6 )




#set up HiC with given AGP
#just because it's stupid to get this from here when it's no doubt been messed up in editing. So delete and reload it directly from the time the SSs were made
agp_given <- copy(orig_agp)
#agp_chr
agp_given %>% setnames("orig_agp_chr","agp_chr")
#add scaffold_length
agp_given <- fai[agp_given,on="scaffold"]
#gaps
setorder(agp_given,"agp_chr","agp_start")
agp_given[ , idx:=2*(1:.N)-1 , by=.(agp_chr) ]
agp_given <- agp_given[ , rbind(.SD,data.table(
  scaffold="gap",
  length=200,
  agp_orientation=NA,
  agp_start=0L
)) , by=.(agp_chr,idx) ]


#orientation NA, 1 or -1
agp_given[ , orientation:=swap(agp_orientation,c("+","-"),c(1,-1)) %>% as.integer ]
agp_given[,agp_orientation:=NULL]

#agp_bin (serves as idx now too???)
agp_given[ scaffold != "gap" , hic_bin := 1:.N , by=.(agp_chr) ]
agp_given[ , gap := (scaffold=="gap") ]
  #kill end gap
agp_given <- agp_given[ , .SD[1:(.N-1)] , by=.(agp_chr) ]
agp_given[ , idx:=NULL ]

#fix agp_start and add agp_end
agp_given[, agp_start := c(1,(cumsum(length)[1:(.N-1)]+1) ) , by="agp_chr" ]
# agp_manual_short[, agp_start := c(1,(cumsum(scaffold_length)[1:(.N-1)]+1) ) , by="agp_chr" ]
agp_given[, agp_end := cumsum(length) , by="agp_chr" ]

setnames(agp_given,c("length"),c("scaffold_length"))
agp_given[ , chr := agp_chr ]
agp_given[ , agp_chr := NULL ]
agp_given <- chrNames_ep(T)[,.(chr,agp_chr)][agp_given,on=.(chr)]

agp_given_bed <- agp_given[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr )]

chrlen <- agp_given[scaffold != "gap" , .(
            length=max(agp_end),
            truechr=T
          ) , by=.(chr)]
chrlen <- chrNames_ep(T)[chrlen,on=.(chr)]
chrun_length <- Sm_v1$info[ ! scaffold %in% agp_given$scaffold ]$length %>% sum
chrlen <- rbind(chrlen,data.table(
            alphachr="Un",
            chr=NA_integer_,
            agp_chr="chrUn",
            length=chrun_length,
            truechr=F
          )
)


Sm_v1$hic_maps$orig_agp$agp <- agp_given
Sm_v1$hic_maps$orig_agp$agp_bed <- agp_given_bed
Sm_v1$hic_maps$orig_agp$chrlen <- chrlen


hic_nucfile <- "data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp_split.nuc.txt"
Sm_v1$hic_maps$orig_agp <- add_psmol_fpairs_ep( assembly = Sm_v1 , hic_map = Sm_v1$hic_maps$orig_agp , nucfile = hic_nucfile ) # adds $links and $frags
Sm_v1$hic_maps$orig_agp$hic_1Mb <- bin_hic_step(hic=Sm_v1$hic_maps$orig_agp$links, frags=Sm_v1$hic_maps$orig_agp$frags, binsize=2e5, chrlen=Sm_v1$hic_maps$orig_agp$chrlen, cores=21)
Sm_v1$hic_maps$orig_agp$hic_1Mb$norm <- normalize_cis(Sm_v1$hic_maps$orig_agp$hic_1Mb, ncores=21, percentile=2 )
#debugonce(find_inversions_ep)
Sm_v1$hic_maps$orig_agp$hic_1Mb$asymmetry <- find_inversions_ep( hic_map=Sm_v1$hic_maps$orig_agp,links=Sm_v1$hic_maps$orig_agp$hic_1Mb$norm , cores=30 )

p0 <- plot_contact_matrices_ep(hicmap=Sm_v1$hic_maps$orig_agp)
p0 <- plot_asymmetry_agp( agp=Sm_v1$hic_maps$orig_agp$agp , asymmetry=Sm_v1$hic_maps$orig_agp$hic_1Mb$asymmetry )
p0
#ggsave(filename = file.path("plots","hic_contact_plot_agp0.pdf") , plot = p0 , device = "pdf" , units = "in" , width = 6 , height = 40 )
#ggsave(filename = file.path("plots","hic_asym_plot_agp0.pdf") , plot = p0 , device = "pdf" , units = "in" , width = 6 , height = 40 )



#check for chimeras
Sm_v1 <- add_hic_cov(Sm_v1, scaffolds=NULL, binsize=1e3, binsize2=1e6, minNbin=50, innerDist=1e5, cores=1)

longest_scaffolds <- fai[order(-length)][1:100]$scaffold

cov_ls <- Sm_v1$cov[longest_scaffolds , on="scaffold"]

dev.off()
pdf(file = "plots/hic_cov_longest100.pdf",paper = "a4")
par(mfrow=c(4,2))
for (s in longest_scaffolds ){
  ce("Scaffold", s)
  p <- cov_ls[scaffold==s][order(bin)]
  plot( p$bin*1e3 , p$n , type='l', lwd=1, col=1  )
  title(s,cex=.5)
}
dev.off()
#





#get marker positions on original scaffolds
markers_to_illumina_scaffs <- fread("data/Eggplant_illumina_scaffold_map_positions.csv")[,ill_scaffold:=sub("_BLOCK_.*","",marker)][]


#get scaffs to superscaffs conversion file
illumina_to_superscaffs <- fread("data/HybridScaffolding.T11.agp",select=c(1:3,5,6,9),col.names=c(
    "scaffold",
    "ss_start",
    "ss_end",
    "ss_type",
    "ill_scaffold",
    "ss_orientation"
  )
)
illumina_to_superscaffs[,ss_length:=max(ss_end),by=.(scaffold)]
illumina_to_superscaffs[ , ss_type:=NULL][]

markers_to_illumina_scaffs[!ill_scaffold %in% illumina_to_superscaffs$ill_scaffold]


# markers_to_illumina_scaffs$ill_scaffold %prop_in% illumina_to_superscaffs$ill_scaffold
# illumina_to_superscaffs[!ill_scaffold %in% markers_to_illumina_scaffs$ill_scaffold]$ill_scaffold
# markers_to_illumina_scaffs[!ill_scaffold %in% illumina_to_superscaffs$ill_scaffold]$ill_scaffold
# markers_to_illumina_scaffs[!ill_scaffold %in% Sm_v1$info$scaffold]
# markers_to_illumina_scaffs[!ill_scaffold %in% Sm_v1$info$scaffold]



markers_to_illumina_scaffs[!grepl("DRAFT",ill_scaffold)]

to_SMEL_DRAFT <- markers_to_illumina_scaffs[ill_scaffold %in% Sm_v1$info$scaffold]
to_Super <- markers_to_illumina_scaffs[!ill_scaffold %in% Sm_v1$info$scaffold]

#checks, does every boy (marker) win a prize (match a Sm_v1 scaffold)? (should all be empty)
to_SMEL_DRAFT[!ill_scaffold %in% Sm_v1$info$scaffold]
to_Super[!ill_scaffold %in% illumina_to_superscaffs$ill_scaffold]

#Create joiner columns
to_SMEL_DRAFT[,scaffold := ill_scaffold]
to_SMEL_DRAFT[ , ss_start := 1L ]
to_SMEL_DRAFT[ , ss_end := NA_integer_ ] # remember to replace with real ends
marker_map <- rbind(to_SMEL_DRAFT,illumina_to_superscaffs[ , .(scaffold,ss_start,ss_end,ill_scaffold) ][to_Super,on=.(ill_scaffold)])

#see what's unjoined. These all match split superscaffolds with a "_a" or such appended
marker_map[!scaffold %in% Sm_v1$info$scaffold]

#find the necessary info for the split ones
add_splitters <- ldply( 1:nrow(marker_map[!scaffold %in% Sm_v1$info$scaffold]) , function(i) { 
  place_me <- marker_map[!scaffold %in% Sm_v1$info$scaffold][i]
  place_me[,frag_end := ss_end]
  setnames(place_me,"scaffold","old_scaffold")
  match_me <- paste0(place_me$old_scaffold,"_.$")
  segments <- Sm_v1$info[grepl(match_me,scaffold)][order(scaffold)]
  segments[,frag_end := cumsum(orig_end)]
  #join placer to end
  segments[place_me,on=.(frag_end),roll=-Inf][,colnames(marker_map),with=F]
})

marker_map <- rbind(marker_map[scaffold %in% Sm_v1$info$scaffold],add_splitters)

#check we got them all now
marker_map[!scaffold %in% Sm_v1$info$scaffold]

#for non-supers, put ss_end as the length./ So that the marker, whose position will be given as the middle, does that correctly.
fai <- fread("data/scaffolds/Scaffold_to_HIC.fa.fai",col.names = c("scaffold","length") , select = 1:2 )
marker_map <- fai[marker_map,on=.(scaffold)]
marker_map[is.na(ss_end),ss_end := length]
marker_map[,pos_in_scaff := pmean(ss_start,ss_end) ]
#marker_map <- marker_map[, chrNames_ep()[,.(gmap_lg=chr,linkage_group=alphachr)][.SD,on="linkage_group"] ]
#add to the list
Sm_v1$marker_map <- marker_map

  #markers <- marker_map
#visualise markers on agp

#debugonce(plot_genetic_map_vs_agp_ep)
plot_genetic_map_vs_agp_ep( markers = Sm_v1$marker_map , agp = Sm_v1$hic_maps$orig_agp$agp )

####################################################################################################################
####################################################################################################################
#saveRDS(Sm_v1,"data/Sm_v1.Rds")
####################################################################################################################
####################################################################################################################


d

p <- d[,.N,by=.(colour,linkage_group)]
plot( x = frank(p$linkage_group,ties.method = "dense") , y = rep(1,nrow(p))  , col = p$colour , cex=4 , pch = 20 )

d[ agp_chr==1 & pos_in_agp < 77e6 ][order(pos_in_agp)]

#Super-Scaffold_719 is notably "out of place" in chr01 ... follow it back
orig_agp <- fread("data/Scaffolds_superscaffolds_eggplant_to_SMELV3.agp",col.names=c("scaffold","agp_orientation","orig_agp_chr","agp_start"))
orig_agp[scaffold=="Super-Scaffold_719"]

illumina_to_superscaffs[scaffold=="Super-Scaffold_719" & ss_orientation != "map"]
illumina_to_superscaffs[scaffold=="Super-Scaffold_719" & ss_orientation != "map"]$ill_scaffold -> scaffs_in_ss719

markers_to_illumina_scaffs[ill_scaffold %in% scaffs_in_ss719]
#ok yeah, it seems to be a genuine problem with this agp (the file or the published version ... ?) YES there was. And maybe still is. And there's essentially fuck-all I can do about it at the moment. Oh well.
d[agp_chr==1 & linkage_group!="E1"]
#
#



#########################################           LOAD FROM HERE        ##########################################
####################################################################################################################
#Sm_v1 <- readRDS("data/Sm_v1.Rds")
#Sm_v2 <- copy(Sm_v1)
####################################################################################################################
####################################################################################################################

#anchor scaffs with popseq. The goal is simply to add a reliable "popseq_cM" column to Sm_v1$info

s <- Sm_v2$marker_map[ !is.na(cM) , if(nu(linkage_group)>1) {
    tbl <- table(linkage_group)
    best <- which(order(-tbl)==1)
    nextbest <- which(order(-tbl)==2)
    .(cM=median(cM),p_lg1_2=tbl[best]/tbl[nextbest],linkage_group=names(tbl)[best],n_markers=.N)
  } else {
    .(cM=median(cM),p_lg1_2=Inf,linkage_group=linkage_group[1],n_markers=.N)
  }
 , by=.(scaffold) ]
s[ p_lg1_2 < 2 , linkage_group := NA ]

Sm_v2$info <- s[ Sm_v2$info , on="scaffold" ]
  Sm_v2$info <- chrNames_gmap_ep()[,.( popseq_chr = chr , linkage_group )][ Sm_v2$info , on="linkage_group" ]
Sm_v2$fragdata$info <- s[Sm_v2$fragdata$info,on="scaffold"]
  Sm_v2$fragdata$info <- chrNames_gmap_ep()[,.( popseq_chr = chr , linkage_group )][ Sm_v2$fragdata$info , on="linkage_group" ]


#suck up a few more scaffolds using Hi-C
Sm_v2 <- anchor_scaffolds_gmap_ep( assembly=Sm_v2 , nlinks_min = 5 , r12_min = 0.6 )
#check
# Sm_v2$fragdata$info[,.N,by=.(scaffold,agp_chr)][,.N,by=.(agp_chr)][order(agp_chr)]
# Sm_v2$info[,.N,by=.(agp_chr)][order(agp_chr)]

#assemble a map for each LG. This is a confused column. Change it here to have the same meaning as it does in the hic_map function (only assigned after hic map construction)
Sm_v2$info[,chr := agp_chr]
Sm_v2$info[, agp_chr := NULL]
Sm_v2$fragdata$info[,chr := agp_chr]
Sm_v2$fragdata$info[,agp_chr := NULL]

chrNames(agp=T)
#remake the hic map
Sm_v2$hic_maps$agp_lgs <- hic_map( info=Sm_v2$fragdata$info , assembly = Sm_v2 , frags=Sm_v2$fragdata$bed )


#########################
saveRDS(Sm_v2,"data/Sm_v2.Rds")
#Sm_v2 <- readRDS("data/Sm_v2.Rds")
###########################

#compare
Sm_v2$hic_maps$agp_lgs$agp[scaffold != "gap"]$scaffold_length %>% sum
Sm_v2$hic_maps$orig_agp$agp[scaffold != "gap"]$scaffold_length %>% sum





plot_genetic_map_vs_agp_ep( markers = Sm_v2$marker_map , agp = Sm_v2$hic_maps$orig_agp$agp )
plot_genetic_map_vs_agp_ep( markers = Sm_v2$marker_map , agp = Sm_v2$hic_maps$agp_lgs$agp  )




hic_nucfile <- "data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp_split.nuc.txt"
Sm_v2$hic_maps$agp_lgs <- add_psmol_fpairs_ep( assembly = Sm_v2 , hic_map = Sm_v2$hic_maps$agp_lgs , nucfile = hic_nucfile ) # adds $links and $frags
Sm_v2$hic_maps$agp_lgs$hic_1Mb <- bin_hic_step(hic=Sm_v2$hic_maps$agp_lgs$links, frags=Sm_v2$hic_maps$agp_lgs$frags, binsize=2e5, chrlen=Sm_v2$hic_maps$agp_lgs$chrlen, cores=21)
Sm_v2$hic_maps$agp_lgs$hic_1Mb$norm <- normalize_cis( Sm_v2$hic_maps$agp_lgs$hic_1Mb, ncores=21, percentile=2 )
Sm_v2$hic_maps$agp_lgs$hic_1Mb$asymmetry <- find_inversions_ep( hic_map=Sm_v2$hic_maps$agp_lgs,links=Sm_v2$hic_maps$agp_lgs$hic_1Mb$norm , cores=30 )

p0 <- plot_contact_matrices_ep(hicmap=Sm_v2$hic_maps$agp_lgs)
p1 <- plot_asymmetry_agp( agp=Sm_v2$hic_maps$agp_lgs$agp , asymmetry=Sm_v2$hic_maps$agp_lgs$hic_1Mb$asymmetry )
ggsave(p0,filename="plots/hic_contact_plot_newanchor_sepLG.pdf",device="pdf",width=5,height=45,units="in")
ggsave(p1,filename="plots/hic_asym_plot_newanchor_sepLG.pdf",device="pdf",width=5,height=25,units="in")


p0 <- plot_contact_matrices_ep(hicmap=Sm_v2$hic_maps$orig_agp)
p1 <- plot_asymmetry_agp( agp=Sm_v2$hic_maps$orig_agp$agp , asymmetry=Sm_v2$hic_maps$orig_agp$hic_1Mb$asymmetry )
ggsave(p1,filename="plots/hic_asym_plot_origAGP.pdf",device="pdf",width=5,height=25,units="in")
ggsave(p0,filename="plots/hic_contact_plot_origAGP.pdf",device="pdf",width=5,height=45,units="in")

#
#
#






#more maps!
maps <- fread("data/maps/allmarkers_to_Scaffold_to_HIC.bl",col.names=c("qseqid", "sseqid" , "slength" , "qlength" , "match_len" , "qstart" , "qend" , "sstart" , "send" , "pct_id" , "evalue" , "bitscore" ))

#Wanna filter hits? Here would be the place.

maps[ , map     := sub("^(.*):::(.*):::(.*):::(.*):::(.*)","\\1",qseqid) ]
maps[ , map_chr := sub("^(.*):::(.*):::(.*):::(.*):::(.*)","\\2",qseqid) ]
maps[ , map_chr := sub(".*(\\d{2})","\\1", map_chr , perl = T ) %>% as.integer  ]
maps[ , cM      := sub("^(.*):::(.*):::(.*):::(.*):::(.*)","\\4",qseqid) %>% as.numeric ]
maps[ , marker      := sub("^(.*):::(.*):::(.*):::(.*):::(.*)","\\5",qseqid)  ]

#now glue them onto the assembly! HHHHERE

#colnames(Sm_v2$marker_map)

Sm_v2$alt_maps <- maps[,.( map , scaffold=sseqid , length=slength , marker , cM , linkage_group=map_chr , pos_in_scaff=pmean(qstart,qend) )] %>%
  rbind( Sm_v2$marker_map[,.( map="BAR2019" , scaffold , length , marker , cM , linkage_group , pos_in_scaff )] )

Sm_v2$alt_maps

markers = Sm_v2$alt_maps
agp = Sm_v2$hic_maps$agp_lgs$agp 

#playground
Sm_v3 <- copy(Sm_v2)

Sm_v3$marker_map <- Sm_v3$alt_maps

Sm_v3$marker_map$cM <- NA

Sm_v3$info$popseq_chr <- NULL
Sm_v3$info$chr <- NULL
Sm_v3$info$best_hic_chr  <- NULL
Sm_v3$info$n_scaff_links <- NULL
Sm_v3$info$scaff_links_r12 <- NULL
Sm_v3$fragdata$info$linkage_group <- NULL
Sm_v3$fragdata$info$chr <- NULL
Sm_v3$fragdata$info$popseq_chr <- NULL

#do initial (naive) anchoring


popseq_assignments <- Sm_v3$marker_map[map != "BAR2019"][  , .( popseq_chr = {
  #if(.GRP==192){browser()}
  tbl <- table(linkage_group)
  tbl <- tbl[rev(order(tbl))]
  lgs <- names(tbl)
  if(length(lgs)==1 | ((tbl[1]/(tbl[1]+tbl[2])) > 0.6)){
    lgs[1] %>% as.integer
  } else {
    NA_integer_
  }
} ) , by=.(scaffold) ]


#join to fraginfo and info
popseq_assignments[Sm_v3$info,on=.(scaffold)] -> Sm_v3$info
popseq_assignments[Sm_v3$fragdata$info,on=.(scaffold)] -> Sm_v3$fragdata$info


Sm_v3 <- anchor_scaffolds_gmap_ep( assembly = copy(Sm_v3) , nlinks_min = 5 , r12_min = 0.6 )


#assemble a map for each LG and fix colnames resulting from crappy programming
Sm_v3$info[,chr := agp_chr]
Sm_v3$info[, agp_chr := NULL]
Sm_v3$fragdata$info[,chr := agp_chr]
Sm_v3$fragdata$info[,agp_chr := NULL]



#remake the hic map
Sm_v3$hic_maps$pub_maps <- hic_map( info=Sm_v3$fragdata$info , assembly = Sm_v3 , frags=Sm_v3$fragdata$bed )

#########################
saveRDS(Sm_v3,"data/Sm_v3.Rds")
#Sm_v3 <- readRDS("data/Sm_v3.Rds")
###########################

plot_genetic_map_vs_agp_ep( markers = Sm_v2$marker_map , agp = Sm_v3$hic_maps$pub_maps$agp )


plot_many_maps_vs_agp_ep( markers = Sm_v3$alt_maps , agp = Sm_v3$hic_maps$pub_maps$agp )
plot_many_maps_vs_agp_ep( markers = Sm_v3$alt_maps , agp = Sm_v2$hic_maps$agp_lgs$agp )



hic_nucfile <- "data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp_split.nuc.txt"
Sm_v3$hic_maps$pub_maps <- add_psmol_fpairs_ep( assembly = Sm_v3 , hic_map = Sm_v3$hic_maps$pub_maps , nucfile = hic_nucfile ) # adds $links and $frags
Sm_v3$hic_maps$pub_maps$hic_1Mb <- bin_hic_step(hic=Sm_v3$hic_maps$pub_maps$links, frags=Sm_v3$hic_maps$pub_maps$frags, binsize=2e5, chrlen=Sm_v3$hic_maps$pub_maps$chrlen, cores=21)
Sm_v3$hic_maps$pub_maps$hic_1Mb$norm <- normalize_cis( Sm_v3$hic_maps$pub_maps$hic_1Mb, ncores=21, percentile=2 )
Sm_v3$hic_maps$pub_maps$hic_1Mb$asymmetry <- find_inversions_ep( hic_map=Sm_v3$hic_maps$pub_maps,links=Sm_v3$hic_maps$pub_maps$hic_1Mb$norm , cores=30 )

p0 <- plot_contact_matrices_ep( hicmap = Sm_v3$hic_maps$pub_maps )
ggsave(p0,filename="plots/hic_contact_plot_pubmap_anchor.pdf",device="pdf",width=5,height=35,units="in")

p1 <- plot_asymmetry_agp( agp=Sm_v3$hic_maps$pub_maps$agp , asymmetry=Sm_v3$hic_maps$pub_maps$hic_1Mb$asymmetry )
ggsave(p1,filename="plots/hic_asym_plot_pubmap_anchor.pdf.pdf",device="pdf",width=5,height=25,units="in")


Sm_v3$hic_maps$pub_maps












#check chr assignments after map-usage (bar19 vs all public maps), and ALSO after hic-curating between these two
popseq_assignments_1 <- Sm_v3$marker_map[map != "BAR2019"][  , .( popseq_chr_pubmaps = most_frequent_by_margin(linkage_group) ) , by=.(scaffold) ]
popseq_assignments_2 <- Sm_v3$marker_map[map == "BAR2019"][  , .( popseq_chr_bar19 = most_frequent_by_margin(linkage_group) ) , by=.(scaffold) ]
cmp_assignments <- popseq_assignments_1[popseq_assignments_2,on="scaffold"]
setorder(cmp_assignments,popseq_chr_pubmaps,popseq_chr_bar19)
setorder(cmp_assignments,popseq_chr_bar19,popseq_chr_pubmaps)
cmp_assignments[,.N,by=.(popseq_chr_pubmaps,popseq_chr_bar19)][!is.na(popseq_chr_pubmaps)] #THIS is the best basis for assignment I say

###HHHHHERE!
#we're going with the bar19 map anchoring, and now to get it down to LGs representing actual chrs. We'll have to wipe out cM dat from the mini-LGs
Sm_v4 <- readRDS("data/Sm_v3.Rds")

Sm_v4$marker_map <- Sm_v4$alt_maps[map=="BAR2019"]

match_lgs <- c("E7","E3","E1","E2A","E11","E6","E12","E10A","E4","E8B","E5","E10B","E9","E2B","E8A")
match_chrs <- c(7,   3,   1,    2,   11,    6,  12,    10,    4,   8,   5  ,  11,   9,    2,    8)
Sm_v4$marker_map[,numeric_chr := swap(linkage_group , match_lgs , match_chrs )  %>% as.integer ]


bar19_scaff_chr_cM <- Sm_v4$marker_map[, .(assigned_chr=most_frequent_by_margin( numeric_chr )) ,by="scaffold"]
t <- bar19_scaff_chr_cM[Sm_v4$marker_map, on="scaffold"]
bar19_scaff_chr_cM <- t[assigned_chr==numeric_chr,.(cM=median(cM)),by=.(scaffold,chr=assigned_chr)]
rm(t)
#
Sm_v4$info[,chr := NULL]
Sm_v4$info[,cM := NULL]
Sm_v4$info <- bar19_scaff_chr_cM[Sm_v4$info,on="scaffold"]
Sm_v4$fragdata$info[,chr := NULL]
Sm_v4$fragdata$info[,cM := NULL]
Sm_v4$fragdata$info <- bar19_scaff_chr_cM[Sm_v4$fragdata$info,on="scaffold"]
Sm_v4$fragdata$info[!is.na(cM)]
Sm_v4$info[!is.na(cM),.N,by=chr]

Sm_v4 <- anchor_scaffolds_gmap_ep( assembly=Sm_v4 , nlinks_min = 5 , r12_min = 0.6 )

Sm_v4$fragdata$info[,chr:=agp_chr]
Sm_v4$fragdata$info[,agp_chr:=NULL]
Sm_v4$info[,chr:=agp_chr]
Sm_v4$info[,agp_chr:=NULL]



#remake the hic map
Sm_v4$hic_maps$bar19_gmap <- hic_map( info=Sm_v4$fragdata$info , assembly = Sm_v4 , frags=Sm_v4$fragdata$bed )
Sm_v4$hic_maps$bar19_gmap$agp[,nu(agp_chr)]
hic_nucfile <- "data/scaffolds/Scaffold_to_HIC_DpnII_fragments_100bp_split.nuc.txt"
Sm_v4$hic_maps$bar19_gmap <- add_psmol_fpairs_ep( assembly = Sm_v4 , hic_map = Sm_v4$hic_maps$bar19_gmap , nucfile = hic_nucfile ) # adds $links and $frags
Sm_v4$hic_maps$bar19_gmap$hic_1Mb <- bin_hic_step(hic=Sm_v4$hic_maps$bar19_gmap$links, frags=Sm_v4$hic_maps$bar19_gmap$frags, binsize=2e5, chrlen=Sm_v4$hic_maps$bar19_gmap$chrlen, cores=21)
Sm_v4$hic_maps$bar19_gmap$hic_1Mb$norm <- normalize_cis( Sm_v4$hic_maps$bar19_gmap$hic_1Mb, ncores=21, percentile=2 )
Sm_v4$hic_maps$bar19_gmap$hic_1Mb$asymmetry <- find_inversions_ep( hic_map=Sm_v4$hic_maps$bar19_gmap,links=Sm_v4$hic_maps$bar19_gmap$hic_1Mb$norm , cores=30 )



Sm_v4$hic_maps$bar19_gmap$agp <- bar19_scaff_chr_cM[,.(scaffold,cM)][Sm_v4$hic_maps$bar19_gmap$agp,on="scaffold"]
Sm_v4$hic_maps$bar19_gmap$agp[,linkage_group:=agp_chr]


#########################
saveRDS(Sm_v4,"data/Sm_v4.Rds")
#Sm_v4 <- readRDS("data/Sm_v4.Rds")
###########################

Sm_v4$hic_maps$agp_lgs <- NULL
Sm_v4$hic_maps$pub_maps <- NULL

Sm_v5 <- copy(Sm_v4)

plot_many_maps_vs_agp_ep( markers = Sm_v5$alt_maps , agp = Sm_v4$hic_maps$bar19_gmap$agp )
#14 * 7 in
plot_contact_matrices_ep( Sm_v5$hic_maps$bar19_gmap )
#90 * 12 in
plot_asymmetry_agp( agp=Sm_v5$hic_maps$bar19_gmap$agp , asymmetry=Sm_v5$hic_maps$bar19_gmap$hic_1Mb$asymmetry  )
#20 * 12 in


make_agp_file( agp=Sm_v5$hic_maps$bar19_gmap$agp,file="plots/manual/agp_out.csv" )
make_agp_file( agp=Sm_v5$hic_maps$bar19_gmap$agp,file="plots/manual/agp_manual.csv" )

Sm_v5$hic_maps$manual <- update_agp_from_file( assembly = Sm_v5 , agp_file = "plots/manual/agp_manual5.csv" )

plot_many_maps_vs_agp_ep( markers = Sm_v5$alt_maps , agp = Sm_v5$hic_maps$manual$agp )
#14 * 7 in
plot_contact_matrices_ep( Sm_v5$hic_maps$manual )
#12 * 90 in
plot_asymmetry_agp( agp=Sm_v5$hic_maps$manual$agp , asymmetry=Sm_v5$hic_maps$manual$hic_1Mb$asymmetry  )
#20 * 12 in
plot_cmp_agps(agp_list = list(Sm_v5$hic_maps$manual$agp,Sm_v5$hic_maps$bar19_gmap$agp,Sm_v5$hic_maps$orig_agp$agp))







#Export and validate
#make fasta file with gap
Sm_v5$hic_maps$manual$agp[scaffold == "gap"][1]$scaffold_length -> gap_length
system("cp  data/scaffolds/Scaffold_to_HIC.fa data/scaffolds/Scaffold_to_HIC_plusgap.fa")
Ns <- rep("N",gap_length) %>% paste0(collapse = "")
system(paste0("printf \">gap\n",Ns,"\" >> data/scaffolds/Scaffold_to_HIC_plusgap.fa"))
system("tail data/scaffolds/Scaffold_to_HIC_plusgap.fa")


#prep bed
agp_bed <- agp[,.(agp_start,agp_end,hic_bin,orientation,scaffold)][Sm_v5$hic_maps$manual$agp_bed,on="scaffold"][,.(scaffold,bed_start,bed_end,score,strand,name=paste0(chr,":",agp_start,":",agp_end,":",orientation,":",hic_bin))]
write.table(agp_bed,file="data/assembly/Sm_v5_agp.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
#run me!
#bedtools getfasta -s -name -fi data/scaffolds/Scaffold_to_HIC_plusgap.fa -bed data/assembly/Sm_v5_agp.bed | perl -e '$lastchr="";$first=1; while(<STDIN>){if($_=~/>1::.*\(([^:]*).*\)/){$chr=$1 ; if($chr ne $lastchr){ if(!$first){print "\n"} print ">",$chr,"\n"; $lastchr=$chr; $first=0} } else {chomp; print} } ' > data/assembly/Sm_pseudomolecules_MTRWv5.fasta


#checks
agp_bed[,sum(bed_end)]+nu(agp$agp_chr) #add a newline at the end of each chromosome
#check equal
#cat data/assembly/Sm_pseudomolecules_MTRWv5.fasta | grep -v '^>' | wc -c


#reconstruct scaffs
agp <- Sm_v5$hic_maps$manual$agp[scaffold != "gap",.(scaffold,)]
fai <- fread("data/scaffolds/Scaffold_to_HIC.fa.fai",col.names = c("scaffold","length") , select = 1:2 )
fai[,idx:=1:.N]

scaffs_from_agp <- agp[fai,on="scaffold"][,.(chr=gsub("chr.*_(.*)","\\1",agp_chr)%>%as.integer,bed_start=agp_start-1,bed_end=agp_end,score=1,strand=swap(orientation,c(1,-1),c("+","-")),name=scaffold)]
write.table(scaffs_from_agp,file="data/assembly/Scaffs_from_agp.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

#co-ersion NAs are good, names chrUn correctly
#bedtools getfasta -s -name -fi data/assembly/Sm_pseudomolecules_MTRWv5.fasta -bed data/assembly/Scaffs_from_agp.bed -name | perl -lnpe 's/1::.*\((.*)\)$/$1/' > TEST_reconstructed_scaffs.fasta
#for i in data/scaffolds/Scaffold_to_HIC.fa TEST_reconstructed_scaffs.fasta; do cat $i | fold -60 > t_$(basename $i); done
#diff t_Scaffold_to_HIC.fa t_TEST_reconstructed_scaffs.fasta -s
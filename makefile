#folder structure must be predone
wd?=/filer-dg/agruppen/seq_shared/MTRW_eggplant_assembly
cutadapt?=/opt/Bio/cutadapt/1.15/bin/cutadapt
minimap?=/opt/Bio/minimap2/2.14/bin/minimap2
samtools?=/opt/Bio/samtools/1.9/bin/samtools
bedtools?=/opt/Bio/bedtools/2.26.0/bin/bedtools
novosort?=/opt/Bio/novocraft/V3.06.05/bin/novosort
bgzip?=/opt/Bio/bcftools/1.9/bin/bgzip

refname?=Scaffold_to_HIC.fa

all :

listify_vcf_for_LD:
	zcat data/map_vcf/All_SNP_HI-C_DP15_missing20.poli.recode_SORTED_maf0.01.vcf.gz | scripts/vcf_listify_for_R.awk > data/map_vcf/All_SNP_HI-C_DP15_missing20.poli.recode_SORTED_maf0.01.tsv

link_hic_data :
	cd $(wd)/data/hic; \
	ln -s /hsm/novaseq/GGR/NOVASEQ/2019/190705_A00550_0039_BHCM2HDRXX/stein/2037510/2/* .

make_scaff_faidx :
	echo "This is done automatically by the digest_ref step, providing bedtools works as usual."

digest_ref :
	cd $(wd)/data/scaffolds; \
	$(wd)/scripts/digest_emboss.zsh --ref $(wd)/data/scaffolds/$(refname) --enzyme 'DpnII' --sitelen 4 --minlen 100 --restrict '/opt/Bio/EMBOSS/6.6.0/bin/restrict' --rebase $(wd)/scripts/link_emboss_e.txt --bedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools' --name MTRW_eggplant_hic; \
	mv *.err $(wd)/err/

minimap_idx :
	$(minimap) -x sr -d $(wd)/data/scaffolds/$(basename $(refname)).mni $(wd)/data/scaffolds/$(refname) 2> $(wd)/err/minimap_indexing.err

map_hic :
	$(wd)/scripts/run_hic_mapping_eggplant.zsh \
	--linker "GATCGATC" \
	--threads 15 \
	--data_dir $(wd)/data/hic \
	--err_dir $(wd)/err \
	--out_dir $(wd)/mapping \
	--tmp $(wd)/temp \
	--name MTRW_eggplant_hic \
	--ref $(wd)/data/scaffolds/$(basename $(refname)).mni \
	--bed $(wd)/data/scaffolds/$(basename $(refname))_DpnII_fragments_100bp.bed \
	--minq 10 \
	--minlen_read 30 \
	--maxlen 500 \
	--onlyfrag 0 \
	--batchmem 5G \
	--indexsize 50G \
	--cutadapt $(cutadapt) \
	--minimap $(minimap) \
	--samtools $(samtools) \
	--bedtools $(bedtools) \
	--novosort $(novosort) \
	--bgzip $(bgzip)
	\
	mv *.err $(wd)/err
	mv MTRW_eggplant_hic* $(wd)/mapping/

check_5p_linker_freq :
	zcat data/hic/2037510_Eggplant_HiC165_S66_L002_R2_001.fastq.gz | head -100000 | awk '(NR%4)==2 {++a[substr($1,0,4)]} END { OFS="\t" ; for (key in a) { print key , a[key] } }' | sort -k2 -n

hic_stats :
	scripts/mapping_stats_ep.zsh --mapdir ./mapping --errdir ./err --prefix MTRW_eggplant_hic > plots/mapping_stats.txt

blast_markers :
	./scripts/blastn_basic.zsh data/scaffolds/Scaffold_to_HIC.fa data/scaffolds/Scaffold_illumina_w_genetmap.fa 5; \
	mv *.bl *_blastn_basic.err mapping



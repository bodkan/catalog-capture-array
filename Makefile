# Every SNC from the catalog of fixed modern human-derived changes will be
# overlapped by two 52bp probes, one will carry a derived allele and second
# will carry an ancestral (archaic-like) allele. Two more probes will be
# flanking the site without any overlap.
#
#                              SNC
#                     25bp      |      26bp
#              <----------------D----------------->
#              <----------------A----------------->
# <----------------------------> <---------------------------->
#              52bp                           52bp
probe_length := 52
 
# input/output directories
scripts_dir := ./scripts
output_dir := ./output
tmp_dir := ./tmp
DIRS := $(output_dir) $(tmp_dir)

script := $(scripts_dir)/calc_probe_coords.py

# input data
catalog_file := /mnt/454/Altaiensis/users/fernando/Combined_Catalog_newf/HumDerived_bothgq30/Genome_VEP.tsv
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa
chrom_info := $(tmp_dir)/chrom_info.txt

# intermediate data
probe_coordinates := $(tmp_dir)/probe_coordinates.bed
human_spec_sites := $(tmp_dir)/snc_positions.bed.gz
overlapping_probes := $(tmp_dir)/overlapping_probes.bed
flanking_probes := $(tmp_dir)/flanking_probes.bed

# final probe sequences
probe_sequences := $(output_dir)/probe_sequences.txt.gz

.PHONY: probes clean

probes: $(DIRS) $(probe_sequences)
	
$(probe_sequences): $(probe_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_tmp -tab
	sed 's/$$/CACTGCGG/' $@_tmp | gzip > $@
	rm $@_tmp

$(probe_coordinates): $(overlapping_probes) $(flanking_probes)
	cat $(overlapping_probes) $(flanking_probes) > $@_tmp
	sort -k1,1n -k2,2n $@_tmp > $@
	rm $@_tmp

$(overlapping_probes): $(human_spec_sites)
	python3 $(script) \
	    --in_file=$< \
	    --out_file=$@ \
	    --probe_length=$(probe_length) \
	    --tiling_step=100 \
	    --flank_length=26

$(flanking_probes): $(human_spec_sites) $(chrom_info)
	bedtools flank -i $(human_spec_sites) -g $(chrom_info) -b $(probe_length) > $@

$(human_spec_sites):
	tail -n+2 $(catalog_file) | \
	grep -v "^X" | \
	awk 'BEGIN {FS = "\t"; OFS = "\t"} \
	     { if (($$21 == "A/A,A/A") && ($$23 > 0.9999)) print $$2, $$3 - 1, $$3 }' | \
	uniq | \
	gzip > $@

# download table of chromosome lengths from UCSC
$(chrom_info):
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz | \
	gunzip | \
	cut -f1,2 | \
	grep -w "chr[X,Y,0-9]*" | \
	sed 's/^chr//' | \
	sort -k1,1V > $@

$(DIRS):
	mkdir $@

clean:
	rm -rf $(DIRS)

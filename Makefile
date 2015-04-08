# default probe-design parameters (in case they are not supplied by the user)
#
#                                  SNC 
#                34bp               |      17bp
# <---------------------------------+---------------->
#                  <----------------+--------------------------------->
#                         17bp      |               34bp
#
probe_length := 52
tiling_step := 17
flank_length := 34

# input/output directories
scripts_dir := ./scripts
output_dir := ./output
tmp_dir := ./tmp
DIRS := $(output_dir) $(tmp_dir)

script := $(scripts_dir)/calc_probe_coords.py

# input data
catalog_file := /mnt/454/Altaiensis/users/fernando/Combined_Catalog_newf/HumDerived_bothgq30/Genome_VEP.tsv
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

# intermediate data
human_spec_sites := $(tmp_dir)/snc_positions.bed.gz
probe_coordinates := $(tmp_dir)/probe_coordinates.bed.gz

# final probe sequences
probe_sequences := $(output_dir)/probe_sequences.txt.gz

.PHONY: probes clean

probes: $(DIRS) $(probe_sequences)
	
$(probe_sequences): $(probe_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_tmp -tab
	sed 's/$$/CACTGCGG/' $@_tmp | tr ':-' '_' | gzip > $@
	rm $@_tmp

$(probe_coordinates): $(human_spec_sites)
	python3 $(script) \
	    --in_file=$(human_spec_sites) \
	    --out_file=$@_tmp \
	    --probe_length=$(probe_length) \
	    --tiling_step=$(tiling_step) \
	    --flank_length=$(flank_length)
	sort -k1,1n -k2,2n $@_tmp | gzip > $@
	rm $@_tmp

$(human_spec_sites):
	tail -n+2 $(catalog_file) | \
	awk 'BEGIN {FS = "\t"; OFS = "\t"} \
	     { if (($$21 == "A/A,A/A") && ($$23 > 0.9999)) print $$2, $$3 - 1, $$3 }' | \
	uniq | \
	gzip > $@

$(DIRS):
	mkdir $@

clean:
	rm -rf $(DIRS)

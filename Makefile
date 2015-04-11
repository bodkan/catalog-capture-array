# Every SNC from the catalog of fixed modern human-derived changes will be
# overlapped by two 52bp probes, one will carry a derived allele and second
# will carry an ancestral (archaic-like) allele. Two more probes will be
# flanking the site without any overlap.
#
#                              SNC
#                     25bp      |      26bp
#              <-----------------D---------------->
#              <--------------------------------->
# <----------------------------> <---------------------------->
#              52bp                           52bp
probe_length := 52
 
# directories
scripts_dir := ./scripts
output_dir := ./output
tmp_dir := ./tmp
DIRS := $(output_dir) $(tmp_dir)

# input data
catalog_file := /mnt/454/Altaiensis/users/fernando/Combined_Catalog_newf/HumDerived_bothgq30/Genome_VEP.tsv
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa
chrom_info := $(tmp_dir)/chrom_info.txt

# intermediate data
human_spec_sites := $(tmp_dir)/snc_positions.bed
overlapping_coordinates := $(tmp_dir)/overlapping_coordinates.bed
flanking_coordinates := $(tmp_dir)/flanking_coordinates.bed
ancestral_sequences := $(tmp_dir)/seq_ancestral_probes.txt
derived_sequences := $(tmp_dir)/seq_derived_probes.txt
flanking_sequences := $(tmp_dir)/seq_flanking_probes.txt

# final probe sequences
final_sequences := $(output_dir)/final_sequences.txt

.PHONY: probes clean

probes: $(DIRS) $(final_sequences)
	
# merge and sort sequences of all generated probe sets and add terminal
# linker sequence to each probe
$(final_sequences): $(derived_sequences) $(ancestral_sequences) $(flanking_sequences)
	cat $^ | \
	    sort -k1,1V -k2,2n | \
	    sed 's/$$/CACTGCGG/' > $@

# get sequences of overlapping probes (the ones that carry fixed modern
# human-specific variants)
$(derived_sequences): $(overlapping_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@ -tab

# get sequences of overlapping probes (the ones that carry ancestral variants
# of SNCs from the catalog)
$(ancestral_sequences): $(overlapping_coordinates)
	# extract sequences of probes overlapping positions from the catalog
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_derived -tab
	
	# get ancestral states at positions of fixed modern human-derived SNCs,
	# add them as another column to the table of sequences above and then
	# substitute the derived state for the ancestral state in each sequence
	# (at position 27, which is the position of SNC in each probe)
	cut -f5 $(human_spec_sites) > ancestral_states.txt
	paste $@_derived ancestral_states.txt | \
		awk -vOFS="\t" '{ \
		    $$2 = substr($$2, 1, 26) $$3 substr($$2, 28); \
	            print $$1, $$2}' \
	        > $@
	
	rm $@_derived ancestral_states.txt

# get sequences of flanking probes
$(flanking_sequences): $(flanking_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@ -tab


# coordinates of probes overlapping SNCs from the catalog
# (probes carrying derived and ancestral variants)
$(overlapping_coordinates): $(human_spec_sites) $(chrom_info)
	bedtools slop -i $(human_spec_sites) -g $(chrom_info) -l 26 -r 25 | \
	    cut -f1,2,3 > $@

# coordinates of probes flanking the SNCs from both sites
$(flanking_coordinates): $(human_spec_sites) $(chrom_info)
	bedtools flank -i $(human_spec_sites) -g $(chrom_info) -b $(probe_length) | \
	    cut -f1,2,3 > $@

# extract genomic coordinates and the derived and ancestral bases for each
# change in the catalog of fixed human-derived changes
$(human_spec_sites):
	tail -n+2 $(catalog_file) | \
	grep -v "^X" | \
	awk 'BEGIN {FS = "\t"; OFS = "\t"} \
	     { if (($$21 == "A/A,A/A") && ($$23 > 0.9999)) print $$2, $$3 - 1, $$3, $$16, $$18 }' | \
	uniq > $@

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

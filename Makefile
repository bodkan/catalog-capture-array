# Every SNC from the catalog of fixed modern human-derived changes will be
# overlapped by two 52bp probes, one will carry a derived allele and second
# will carry an ancestral (archaic-like) allele. Two more probes will be
# flanking the site without any overlap.
#
#                              SNC
#                     26bp      |      25bp
#             <-----------------D---------------->
#             <-----------------A---------------->
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
chr_lengths := $(tmp_dir)/chr_lengths.txt

# coordinates of fixed modern human-derived sites from the catalog
human_derived_sites := $(tmp_dir)/fixed_human_derived_sites.bed

# coordinates of overlapping and flanking probes
overlapping_coordinates := $(tmp_dir)/overlapping_coordinates.bed
flanking_coordinates := $(tmp_dir)/flanking_coordinates.bed

# sequences of overlapping (carrying ancestral or derived variants)
# and flanking probes
ancestral_sequences := $(tmp_dir)/seq_ancestral_probes.txt
derived_sequences := $(tmp_dir)/seq_derived_probes.txt
flanking_sequences := $(tmp_dir)/seq_flanking_probes.txt

# sequences of the final probe set
final_sequences := $(output_dir)/final_sequences.txt

.PHONY: probes clean

probes: $(DIRS) $(final_sequences)
	
# merge and sort sequences of all generated probe sets and add a linker
# sequence to the end of each probe
$(final_sequences): $(derived_sequences) $(ancestral_sequences) $(flanking_sequences)
	cat $^ | \
	    sort -k1,1V -k2,2n | \
	    sed 's/$$/CACTGCGG/' | \
	    uniq > $@

# get sequences of overlapping probes carrying fixed human-derived variants
$(derived_sequences): $(overlapping_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_tmp -tab
	sed 's/\t/d\t/g' $@_tmp > $@
	rm $@_tmp

# get sequences of overlapping probes carrying ancestral variants of SNCs
# that are fixed derived in present-day humans
$(ancestral_sequences): $(overlapping_coordinates)
	# extract sequences of probes overlapping positions from the catalog
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_derived -tab
	
	# get ancestral states at positions of fixed derived SNCs...
	cut -f5 $(human_derived_sites) > ancestral_states.txt
	# ... add them as another column to the table of sequences and
	# substitute the derived state for the ancestral state in each sequence
	# in the second column (at position 27, which is the position of SNC
	# in each probe)
	paste $@_derived ancestral_states.txt | \
		awk -vOFS="\t" '{ \
		    $$2 = substr($$2, 1, 26) $$3 substr($$2, 28); \
	            print $$1"a", $$2}' \
	        > $@
	
	rm $@_derived ancestral_states.txt

# get sequences of flanking probes
$(flanking_sequences): $(flanking_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@ -tab


# coordinates of probes overlapping SNCs from the catalog
# (probes carrying derived and ancestral variants)
$(overlapping_coordinates): $(human_derived_sites) $(chr_lengths)
	bedtools slop -i $(human_derived_sites) -g $(chr_lengths) -l 26 -r 25 | \
	    cut -f1,2,3 > $@

# coordinates of probes flanking the SNCs from both sites
$(flanking_coordinates): $(human_derived_sites) $(chr_lengths)
	bedtools flank -i $(human_derived_sites) -g $(chr_lengths) -b $(probe_length) | \
	    cut -f1,2,3 > $@

# extract genomic coordinates and the derived and ancestral bases for each
# change in the catalog of fixed human-derived changes
$(human_derived_sites):
	tail -n+2 $(catalog_file) | \
	grep -v "^X" | \
	awk 'BEGIN {FS = "\t"; OFS = "\t"} \
	     { if (($$21 == "A/A,A/A") && ($$23 > 0.9999)) print $$2, $$3 - 1, $$3, $$16, $$18 }' | \
	uniq > $@

# download table of chromosome lengths from UCSC
$(chr_lengths):
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

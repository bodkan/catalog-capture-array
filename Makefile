# Every SNP from the specified BED file will be overlapped by two 52bp probes,
# one will carry a reference allele and the second will carry an alternative
# allele. Two more probes will be flanking the site without any overlap.
#
#                              SNC
#                     26bp      |      25bp
#             <-----------------R---------------->
#             <-----------------A---------------->
# <----------------------------> <---------------------------->
#              52bp                           52bp

probe_length := 52
 
# directories
scripts_dir := ./scripts
output_dir := ./output
tmp_dir := ./tmp
DIRS := $(output_dir) $(tmp_dir)

# lengths of chromosomes required for some bedtools utils
chr_lengths := $(tmp_dir)/chr_lengths.txt

# coordinates of overlapping and flanking probes
overlapping_coordinates := $(tmp_dir)/overlapping_coordinates.bed
flanking_coordinates := $(tmp_dir)/flanking_coordinates.bed

# sequences of overlapping and flanking probes
reference_sequences := $(tmp_dir)/seq_reference_probes.txt
alternative_sequences := $(tmp_dir)/seq_alternative_probes.txt
flanking_sequences := $(tmp_dir)/seq_flanking_probes.txt

# sequences of the final probe set
final_sequences := $(output_dir)/final_sequences.txt

.PHONY: probes clean

default:
	@echo "Usage:\n\tmake probes snp_positions=<input BED> ref_genome=<FASTA reference genome>\n"

probes: $(DIRS) $(final_sequences)
	
# merge and sort sequences of all generated probe sets and add a linker
# sequence to the end of each probe
$(final_sequences): $(reference_sequences) $(alternative_sequences) $(flanking_sequences)
	cat $^ | \
	    sort -k1,1V -k2,2n | \
	    sed 's/$$/CACTGCGG/' | \
	    uniq > $@

# get sequences of overlapping probes carrying the reference variants
$(reference_sequences): $(overlapping_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_tmp -tab
	sed 's/\t/d\t/g' $@_tmp > $@
	rm $@_tmp

# get sequences of overlapping probes carrying the alternative variants
$(alternative_sequences): $(overlapping_coordinates)
	# extract sequences of probes overlapping positions from the catalog
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_tmp -tab
	
	# get the alternative states at positions of interest...
	cut -f5 $(snp_positions) > alternative_states.txt
	# ... add them as another column to the table of sequences and
	# substitute the reference allele for the alternative in each sequence
	# in the second column (at position 27, which is the position of SNP
	# in each probe)
	paste $@_tmp alternative_states.txt | \
		awk -vOFS="\t" '{ \
		    $$2 = substr($$2, 1, 26) $$3 substr($$2, 28); \
	            print $$1"a", $$2}' \
	        > $@
	
	rm $@_tmp alternative_states.txt

# get sequences of flanking probes
$(flanking_sequences): $(flanking_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@ -tab


# coordinates of probes overlapping the SNPs of interest
# (probes carrying reference and alternative states)
$(overlapping_coordinates): $(snp_positions) $(chr_lengths)
	bedtools slop -i $(snp_positions) -g $(chr_lengths) -l 26 -r 25 | \
	    cut -f1,2,3 > $@

# coordinates of probes flanking the SNPs from both sides
$(flanking_coordinates): $(snp_positions) $(chr_lengths)
	bedtools flank -i $(snp_positions) -g $(chr_lengths) -b $(probe_length) | \
	    cut -f1,2,3 > $@

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

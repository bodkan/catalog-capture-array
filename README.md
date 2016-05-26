## Design probes to capture a given set of SNPs

Every SNP from a specified input BED file will be overlapped by two 52bp probes,
one carrying a reference allele (R) and the second carrying an alternative allele (A).
Two more probes will be flanking the site without any overlap.

```
                              SNC
                     26bp      |      25bp
             <-----------------R---------------->
             <-----------------A---------------->
 <----------------------------> <---------------------------->
              52bp                           52bp
```

The input has to be in a BED format with five columns, three specifying the coordinates
as usual (chromosome, 0-based position, 1-based position) and the other two columns being
the reference and alternative alleles. Example:

```
1       1136204 1136205 T       C
1       1198046 1198047 C       T
1       1203167 1203168 G       A
1       1501063 1501064 T       G
1       1840017 1840018 C       G
1       1847134 1847135 T       C
1       1865119 1865120 C       G
[...]
```

### How to use it
```
$> make
Usage:
        make probes snp_positions=<input BED> ref_genome=<FASTA reference genome>
```
The input BED file and the reference genome must be specified by assigning to variables `snp_positions`
and `ref_genome`, respectively, while invoking `make`.

### Requirements

* Unix-like system with GNU Make, GNU coreutils, GNU awk and GNU grep (although I don't think there's a reason why this wouldn't work with non-GNU version of these).
* curl
* bedtools (developed and tested on v2.24.0).
* reference genome in the FASTA format

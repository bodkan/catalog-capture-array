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
the reference and alternative alleles.

```
$> make
Usage:
        make probes snp_positions=<path to BED file with SNP positions>
```

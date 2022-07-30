# vcfverifier

Checks that a given VCF file matches a given assembly in FASTA format by
checking that the REF column matches the FASTA file for each record in the
FASTA file (case insensitive)

## Install

First install rust, probably with rustup https://rustup.rs/

Then

```
cargo install --git https://github.com/cmdcolin/vcfverifier
```

## Usage

```
## Generated FASTA index (fai)
samtools faidx myfile.fa

## Run the verifier
vcfverifier --fasta myfile.fa --vcf myfile.vcf.gz
```

Allows plaintext, gzip, or bgzip vcf files as input to the --vcf flag

## Approx speed

Processing chr1 (6.5M rows) of the 1000 genomes dataset takes ~24seconds

```
$ time vcfverifier --fasta hs37d5.fa --vcf ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
Lines processed: 6468347
No mismatching lines found
vcfverifier --fasta ~/Downloads/hs37d5.fa --vcf   24.07s user 0.26s system 99% cpu 24.330 total

```

## Note

My first rust project!

Uses faimm to memory-map the indexed FASTA file, keeping memory usage low (the
entire FASTA does not have to be loaded into memory and the VCF is read line by
line) https://github.com/veldsla/faimm

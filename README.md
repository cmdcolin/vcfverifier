# vcfverifier

Checks that a given VCF file matches a given assembly in FASTA format

## Usage

```
## Generated FASTA index (fai)
samtools faidx myfile.fa

## Run the verifier
vcfverifier --fasta myfile.fa --vcf myfile.fa.gz
```

## Note

My first rust project!

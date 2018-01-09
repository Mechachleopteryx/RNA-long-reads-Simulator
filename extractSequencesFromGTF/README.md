The gffread command line is:
```
gffread/gffread -g <reference.fasta> -w <transcripts.fa> <gtf_file>
```

We just need to take extra care that the chromosomes names in the genome_dir is the same chromosome names in the gtf file

For the sample example, GTFs downloaded from https://www.ensembl.org/info/data/ftp/index.html. Sequences downloaded also from Ensemble.


To run on the sample example (to be used as an example for the real command line):
```
gunzip -k sample_example/Homo_sapiens.GRCh38.91.Y.gtf.gz
gunzip -k sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
gffread/gffread -g sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa -w transcripts.fa sample_example/Homo_sapiens.GRCh38.91.Y.gtf
rm sample_example/Homo_sapiens.GRCh38.91.Y.gtf
rm sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa
rm sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.fai
```

The transcripts will be in file transcripts.fa

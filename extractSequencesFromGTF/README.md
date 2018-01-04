The gffread command line is:
```
gffread/gffread -g <genome_dir> -w <transcripts.fa> <gtf_file>
```

We just need to take extra care that the chromosomes names in the genome_dir is the same chromosome names in the gtf file

For the sample example, GTFs downloaded from https://www.ensembl.org/info/data/ftp/index.html. Sequences downloaded also from Ensemble.


To run on the sample example (to be used as an example for the real command line):
```
gunzip -k sample_example/Homo_sapiens.GRCh38.91.Y.gtf.gz
gunzip -k sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
mv sample_example/Homo_sapiens.GRCh38.dna.chromosome.Y.fa sample_example/Y.fa
gffread/gffread -g sample_example/ -w transcripts.fa sample_example/Homo_sapiens.GRCh38.91.Y.gtf
rm sample_example/Homo_sapiens.GRCh38.91.Y.gtf
rm sample_example/Y.fa
```

The transcripts will be in file transcripts.fa

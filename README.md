# RNA ONT MinION simulator

A tool for simulating RNA Oxford Nanopore Technologies long reads.

This tool relies on FluxSimulator (http://sammeth.net/confluence/display/SIM/Home), distributed as binary in the `thirdparty` folder.

## Dependencies:
```
python 2.7
virtualenv
g++
java
```

## Installing:
```
git clone --recursive https://gitlab.inria.fr/marchet/RNA-long-reads-Simulator
cd RNA-long-reads-Simulator
bash install.sh
```

## Testing:
```
cd <path-to-RNA-long-reads-Simulator>
source venv/bin/activate
cd sample_test/
tar -zxvf sample_test.tar.gz
python ../RNACreator.py -g Mus_musculus.GRCm38.91.chr19.gtf -b gmap_CB_1Donly_to_GRCm38_chr19.sorted.bam -i gmap_CB_1Donly_to_GRCm38_chr19.sorted.bam.bai -r Mus_musculus.GRCm38.dna.chromosome.19.fa
deactivate
```

## Parameters:
```
usage: RNACreator.py [-h] [-g GTFFILEPATH] [-b BAMFILEPATH] [-i BAIFILEPATH]
                     [-r GENOMEREFPATH] [-c COVERAGE] [-o OUTPUTDIRPATH]
                     [--version]

RNACreator - Simulation of whole transcriptome data sets of ONT long reads

optional arguments:
  -h, --help        show this help message and exit
  -g GTFFILEPATH    Path to the GTF file
  -b BAMFILEPATH    Path to the BAM file - used only to produced the error
                    profile
  -i BAIFILEPATH    Path to the BAI file - used only to produced the error
                    profile
  -r GENOMEREFPATH  Path to the reference genome file
  -c COVERAGE       An integer that represents the desired coverage
                    (default=1)
  -o OUTPUTDIRPATH  Path to the output directory (default: .)
  --version         show program's version number and exit
```
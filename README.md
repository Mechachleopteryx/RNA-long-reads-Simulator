# RNA ONT MinION simulator

## Installing:
```
git clone --recursive https://gitlab.inria.fr/marchet/RNA-long-reads-Simulator
cd RNA-long-reads-Simulator
bash install.sh
```

## Testing:
```
cd sample_test/
tar -zxvf sample_test.tar.gz
python ../RNACreator.py -g Mus_musculus.GRCm38.91.chr19.gtf -b gmap_CB_1Donly_to_GRCm38_chr19.sorted.bam -i gmap_CB_1Donly_to_GRCm38_chr19.sorted.bam.bai -r Mus_musculus.GRCm38.dna.chromosome.19.fa --fsPath /data2/leandro/flux_simulator/flux-simulator-1.2.1/bin/flux-simulator
```


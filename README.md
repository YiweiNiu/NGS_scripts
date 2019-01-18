# Scripts used for genomics and bioinformatics


## batchPreprocess.HT.py - submit jobs to the TORQUE cluter in batch

Prepare working directory like this:

```bash
WORKDIR
├── sh
│   ├── batchPreprocess.HT.py
├── jobfiles
├── logs
├── rawsra
│   ├── dataset1
│   │   ├── SRR1
│   │   └── SRR2
│   └── dataset2
│       └── SRR3
├── rawfastq
│   ├── dataset1
│   └── dataset2
├── bowtie2
│   ├── dataset1
│   │   ├── SRR1.bam
│   │   └── SRR2.bam
│   └── dataset2
└── STAR
    ├── dataset1
    │   └── SRR1.bam
    └── dataset2
        └── SRR2.bam
├── tmp

```

Usage:

```bash
python batchPreprocess.HT.tmp.py -h
usage: batchPreprocess.HT.tmp.py [-h] [--version] {rnaseq,atacseq,chipseq} ...

batchPreprocess.HT.tmp.py -- submit jobs to the clusters in batch.

positional arguments:
  {rnaseq,atacseq,chipseq}
    atacseq             Main Function: preprocess ATAC-seq data.
    rnaseq              Main Function: preprocess RNA-seq data.
    chipseq             Main Function: preprocess ChIP-seq data.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each command, type: batchPreprocess.HT.py COMMAND -h
```

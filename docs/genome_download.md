# Downloading reference genomes

## GRCh38 (Human, Ensembl 115)

```bash
mkdir -p ref && cd ref
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

## GRCm39 (Mouse, Ensembl 115)

```bash
mkdir -p ref && cd ref
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
```

## Building the kallisto index

```bash
kb ref \
    -i ref/index.idx \
    -g ref/t2g.txt \
    -f1 ref/cdna.fa \
    ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    ref/Homo_sapiens.GRCh38.115.gtf.gz
```

The pipeline can do this for you via `--genome_fasta` + `--genome_gtf`, or use
`--prebuilt_index` + `--t2g` to skip the build.

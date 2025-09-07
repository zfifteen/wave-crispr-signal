# Dataset Provenance for Signal-Theoretic CRISPR Experiment

## Primary Datasets

### Doench 2016 Human CRISPR Efficiency Dataset
- **Name**: doench2016.csv
- **Source**: Doench et al. (2016) "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9"
- **URL**: https://www.nature.com/articles/nbt.3437
- **License**: Available for research use
- **Organism**: Homo sapiens
- **Type**: On-target efficiency regression
- **Samples**: 583 human CRISPR guide sequences
- **SHA256**: 6c2c38934aa771e74ecb42a7b477e0541694aba805ab5de529e87448416ec892
- **Local Path**: data/doench2016.csv

### BioGRID-ORCS Human Dataset  
- **Name**: BioGRID-ORCS Homo sapiens v1.1.17
- **Source**: BioGRID Open Repository of CRISPR Screens
- **URL**: https://orcs.thebiogrid.org/
- **License**: Open data license
- **Citation**: Oughtred et al. (2021) "The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions"
- **Organism**: Homo sapiens
- **Version**: 1.1.17
- **Local Path**: data/BIOGRID-ORCS-ALL-homo-sapiens-1.1.17.screens

## Reference Genomes

### Human Genome Reference (hg38)
- **Name**: GRCh38/hg38
- **Source**: UCSC Genome Browser
- **URL**: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
- **Build**: GRCh38
- **Files**: 
  - hg38.fa.gz (Complete genome sequence)
  - hg38.chrom.sizes (Chromosome sizes)
- **Download Script**: data/get_hg38/get_hg38.sh
- **Documentation**: data/get_hg38/get_hg38_README.md

## Data Acquisition

To download the reference genome data:

```bash
cd data/get_hg38
bash get_hg38.sh
```

This will download and organize the hg38 reference files according to the structure documented in get_hg38_README.md.

## Validation

All datasets undergo validation through the HumanFASTAValidator to ensure:
- Human organism only (Homo sapiens)
- Valid nucleotide alphabet (A/C/G/T/N)
- Proper sequence format and structure
- Ethics compliance (public datasets only)

## Checksums Verification

To verify data integrity:

```bash
# Verify main dataset
echo "6c2c38934aa771e74ecb42a7b477e0541694aba805ab5de529e87448416ec892  data/doench2016.csv" | sha256sum -c

# Generate checksums for all data files
find data/ -type f -name "*.csv" -o -name "*.fasta" -o -name "*.fa" | xargs sha256sum
```
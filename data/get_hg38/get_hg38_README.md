This document outlines the structure and contents of the `data/get_hg38` directory, which is designed to download and organize human genome reference data (hg38) for analysis.

### Data Acquisition Script

*   **File:** `get_hg38.sh`
*   **Purpose:** This shell script automates the download of specified hg38 data files from the UCSC Genome Browser's FTP server.
*   **Functionality:**
    1.  It checks a `zips/` subdirectory for existing downloaded files to avoid re-downloading.
    2.  For each required file, it downloads the compressed archive from UCSC.
    3.  It creates a new directory named after the archive (e.g., `hg38.fa` for `hg38.fa.gz`).
    4.  It extracts the contents of the archive into the corresponding new directory.
    5.  The original downloaded archive is moved into the `zips/` directory for caching and organization.

### Directory and File Contents

The script generates the following directory structure:

*   `get_hg38.sh`: The executable script to download all data.
*   `zips/`: This directory acts as a cache, storing the original compressed files (`.gz`, `.tar.gz`, `.bw`) downloaded from UCSC.
*   `hg38.fa/`: Contains `hg38.fa`, the complete hg38 human genome sequence in a single FASTA file.
*   `hg38.fa.masked/`: Contains `hg38.fa.masked`, a version of the genome where repeats are soft-masked (represented by lowercase letters).
*   `hg38.chromFa/`: Contains individual FASTA files for each chromosome and contig (e.g., `chr1.fa`, `chrX.fa`).
*   `hg38.chromFaMasked/`: Contains individual repeat-masked FASTA files for each chromosome.
*   `mrna.fa/`: Contains `mrna.fa`, a FASTA file of mRNA sequences from GenBank.
*   `upstream1000.fa/`, `upstream2000.fa/`, `upstream5000.fa/`: These directories contain FASTA files with sequences 1000, 2000, or 5000 base pairs upstream of annotated transcription start sites, respectively. They are used for promoter analysis.
*   `hg38.chrom.sizes`: A two-column text file listing each chromosome and its total length in base pairs.
*   `hg38.gc5Base.wigVarStep/`: Contains `hg38.gc5Base.wigVarStep`, a Wiggle format file detailing GC percentage in 5-base windows across the genome. The original download is a more efficient BigWig (`.bw`) file.
*   `hg38.trf.bed/`: Contains `hg38.trf.bed`, a BED file that annotates the locations of tandem repeats.
*   `hg38.fa.out/`: Contains `hg38.fa.out`, the detailed annotation output from the RepeatMasker tool.
*   `hg38.fa.align/`: Contains `hg38.fa.align`, a file showing the alignments of repetitive elements to the consensus repeat sequences.

To use this data, first execute the `get_hg38.sh` script to ensure all necessary files are downloaded and placed in their correct locations.
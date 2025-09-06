#!/bin/bash

# Shell script to download required hg38 files from UCSC bigZips,
# extract contents into subdirectories named after the files (without extensions),
# and move the original downloaded files to a 'zips' folder.
# Based on download instructions from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
# Uses wget with ftp protocol as recommended.
# Handles .tar.gz (extract with tar), .gz (uncompress with gunzip), and plain files (copy).
# Skips refMrna.fa.gz as it is not available in hg38 bigZips (available in hg19).

# Create zips directory if it doesn't exist
mkdir -p zips

# List of required files (confirmed available in hg38 bigZips)
files=(
  "hg38.fa.gz"
  "hg38.fa.masked.gz"
  "hg38.chromFa.tar.gz"
  "hg38.chromFaMasked.tar.gz"
  "mrna.fa.gz"
  "upstream1000.fa.gz"
  "upstream2000.fa.gz"
  "upstream5000.fa.gz"
  "hg38.chrom.sizes"
  "hg38.gc5Base.bw"
  "hg38.gc5Base.wigVarStep.gz"
  "hg38.trf.bed.gz"
  "hg38.fa.out.gz"
  "hg38.fa.align.gz"
)

base_url="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/"

for file in "${files[@]}"; do
  # Download the file using wget with timestamping
  wget --timestamping "${base_url}${file}" -O "${file}"

  # Determine subdirectory name (remove .tar.gz or .gz)
  if [[ $file == *.tar.gz ]]; then
    subdir="${file%.tar.gz}"
  elif [[ $file == *.gz ]]; then
    subdir="${file%.gz}"
  else
    subdir="${file}"
  fi

  # Create the subdirectory
  mkdir -p "${subdir}"

  # Extract or copy contents
  if [[ $file == *.tar.gz ]]; then
    tar -xzf "${file}" -C "${subdir}"
  elif [[ $file == *.gz ]]; then
    gunzip -c "${file}" > "${subdir}/${file%.gz}"
  else
    cp "${file}" "${subdir}/"
  fi

  # Move the original file to zips
  mv "${file}" zips/
done

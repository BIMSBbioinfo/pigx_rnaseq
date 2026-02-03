#!/usr/bin/env bash
set -euo pipefail

# Cleans FASTA headers to retain only Ensembl transcript IDs, ensuring
# consistency between cDNA FASTA files and transcript_id fields in Ensembl GTF files.

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <input.fa|input.fasta>" >&2
  exit 1
fi

in="$1"

case "$in" in
  *.fa)     base="${in%.fa}";     ext="fa" ;;
  *.fasta) base="${in%.fasta}";  ext="fasta" ;;
  *)
    echo "Error: input file must have .fa or .fasta extension" >&2
    exit 1
    ;;
esac

out="${base}.cleaned.${ext}"

sed '/^>/ s/^\(>[^.]*\)\..*/\1/' "$in" > "$out"

echo "New fasta file is printed here ${out}"

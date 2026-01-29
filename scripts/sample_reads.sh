#!/usr/bin/env bash
# Usage:
#   chmod +x sample_reads.sh
#   ./sample_reads.sh <input_folder> <output_folder> <num_lines> <extension>
#
# Example:
#   ./sample_reads.sh raw_fastqs sampled_reads 100000 fastq.gz

in="$1"; out="$2"; n="$3"; ext="$4"
mkdir -p "$out"
find "$in" -name "*.${ext}" -type f | while read -r f; do
  b=$(basename "$f" ".${ext}")
  echo "Extracting ${n} lines from ${f}"
  zcat "$f"  2>/dev/null | head -n "$n" | gzip -c > "$out/${b}.${ext}"
done


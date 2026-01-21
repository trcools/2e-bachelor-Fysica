#!/usr/bin/env bash
set -euo pipefail

# Build PDF from Markdown using Pandoc + tectonic
cd "$(dirname "$0")"

IN="Algemene-Samenvatting-Kwantum.md"
OUT="Algemene-Samenvatting-Kwantum.pdf"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "Error: pandoc is not installed. On macOS: brew install pandoc" 1>&2
  exit 1
fi

if ! command -v tectonic >/dev/null 2>&1; then
  echo "Error: tectonic is not installed. On macOS: brew install tectonic" 1>&2
  exit 1
fi

pandoc "$IN" \
  --from=markdown+tex_math_dollars \
  --pdf-engine=tectonic \
  --toc \
  --number-sections \
  --include-in-header=preamble.tex \
  -o "$OUT"

echo "Built $OUT"
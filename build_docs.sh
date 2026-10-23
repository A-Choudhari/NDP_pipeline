#!/usr/bin/env bash
# Rebuild NDP_instructions.pdf from NDP_instructions.md (pandoc -> HTML -> headless Chrome).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
pandoc NDP_instructions.md -o NDP_instructions.html --standalone --toc --toc-depth=2 \
    --css=docs/print.css --embed-resources
google-chrome --headless --disable-gpu --no-sandbox --no-pdf-header-footer \
    --print-to-pdf=NDP_instructions.pdf --virtual-time-budget=10000 NDP_instructions.html
rm -f NDP_instructions.html
echo "NDP_instructions.pdf  ($(pdfinfo NDP_instructions.pdf | awk '/^Pages/{print $2}') pages)"

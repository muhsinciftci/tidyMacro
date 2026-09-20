#!/usr/bin/env bash
# Post-render hook: mirror the Quarto output into ../docs, which is the
# directory the published website is served from.
#
# Quarto refuses to manage an output-dir that sits outside its own project
# directory, so it renders into _site/_output and this script mirrors that
# tree.  --delete means docs/ never accumulates pages whose source .qmd has
# been renamed or removed; the excludes protect files a host may need that
# Quarto does not generate.
set -euo pipefail

root="${QUARTO_PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
src="$root/_output/"
dst="$root/../docs/"

mkdir -p "$dst"
rsync -a --delete \
      --exclude '.nojekyll' \
      --exclude 'CNAME' \
      --exclude '.git/' \
      "$src" "$dst"

echo "post-render: synced $(find "$src" -name '*.html' | wc -l | tr -d ' ') pages -> docs/"

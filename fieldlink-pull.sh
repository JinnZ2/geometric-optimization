#!/usr/bin/env bash
# fieldlink-pull.sh — Pull Rosetta-Shape-Core data into .fieldlink/merge_stage/
set -euo pipefail

CONFIG=".fieldlink.json"
STAGE_DIR=".fieldlink/merge_stage/rosetta"
REPO_URL="https://github.com/JinnZ2/Rosetta-Shape-Core.git"
TMP_DIR=$(mktemp -d)

trap 'rm -rf "$TMP_DIR"' EXIT

echo "==> Pulling Rosetta-Shape-Core..."
git clone --depth 1 "$REPO_URL" "$TMP_DIR" 2>/dev/null

mkdir -p "$STAGE_DIR/shapes" "$STAGE_DIR/bridges" "$STAGE_DIR/ontology" "$STAGE_DIR/schema"

cp "$TMP_DIR"/shapes/*.json       "$STAGE_DIR/shapes/"       2>/dev/null || true
cp "$TMP_DIR"/bridges/rosetta-bridges.json "$STAGE_DIR/bridges/" 2>/dev/null || true
cp "$TMP_DIR"/ontology/_vocab.json        "$STAGE_DIR/ontology/" 2>/dev/null || true
cp "$TMP_DIR"/ontology/capabilities.json  "$STAGE_DIR/ontology/" 2>/dev/null || true
cp "$TMP_DIR"/schema/shape.schema.json    "$STAGE_DIR/schema/"   2>/dev/null || true

echo "==> Staged to $STAGE_DIR"
ls -R "$STAGE_DIR"
echo "==> Fieldlink pull complete."

#!/usr/bin/env bash
# Install the MDPA extension into VS Code (server or desktop) without vsce.
# Run from any directory: bash scripts/tools/vscode-mdpa/install.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXT_NAME="kratos-multiphysics.vscode-mdpa-0.1.0"

# Find the extensions folder (server > desktop > fallback).
if [[ -d "$HOME/.vscode-server/extensions" ]]; then
  EXT_DIR="$HOME/.vscode-server/extensions/$EXT_NAME"
elif [[ -d "$HOME/.vscode/extensions" ]]; then
  EXT_DIR="$HOME/.vscode/extensions/$EXT_NAME"
else
  echo "Could not find a VS Code extensions folder under ~/.vscode-server or ~/.vscode" >&2
  exit 1
fi

echo "→ Building extension..."
(cd "$SCRIPT_DIR" && npm run compile)

echo "→ Installing to $EXT_DIR"
rm -rf "$EXT_DIR"
mkdir -p "$EXT_DIR"

# Copy the files VS Code needs (no node_modules, no dev sources).
cp "$SCRIPT_DIR/package.json"              "$EXT_DIR/"
cp "$SCRIPT_DIR/language-configuration.json" "$EXT_DIR/"
cp -r "$SCRIPT_DIR/syntaxes"               "$EXT_DIR/"
cp -r "$SCRIPT_DIR/dist"                   "$EXT_DIR/"
cp -r "$SCRIPT_DIR/media"                  "$EXT_DIR/"

echo "✓ Extension installed."
echo ""
echo "  Reload VS Code (Ctrl+Shift+P → Developer: Reload Window), then"
echo "  open any .mdpa file and use the 'Open MDPA Preview' button or command."

#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$PROJECT_DIR/build"
BIN_DIR="$PROJECT_DIR/bin"

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR" "$BIN_DIR"

cd "$BUILD_DIR"
cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="$BIN_DIR" ..
cmake --build . -j

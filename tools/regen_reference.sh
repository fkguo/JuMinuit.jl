#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-2.1-or-later
#
# Regenerate test/reference_data/*.json from the pinned C++ Minuit2.
#
# Idempotent: re-runs build then harness, overwriting any prior JSON.
# Bit-exact reproducibility requires the same machine/BLAS combination
# documented in test/reference_data/_machine.txt (recorded by this
# script).
#
# Usage:
#   tools/regen_reference.sh

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="$REPO_ROOT/tools"
BUILD_DIR="$TOOLS_DIR/build"
REF_DIR="$REPO_ROOT/test/reference_data"
REF_CPP="$REPO_ROOT/reference/Minuit2_cpp"

# Sanity: reference checkout present
if [[ ! -d "$REF_CPP" ]]; then
    echo "ERROR: $REF_CPP missing — re-clone per docs/UPSTREAM.md" >&2
    exit 1
fi

mkdir -p "$BUILD_DIR" "$REF_DIR"

echo ">>> Configuring CMake (this also configures Minuit2 standalone)"
cmake -S "$TOOLS_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release

echo ">>> Building cpp_trace_harness (also builds Minuit2.a)"
cmake --build "$BUILD_DIR" --target cpp_trace_harness --parallel

echo ">>> Running harness"
"$BUILD_DIR/cpp_trace_harness" "$REF_DIR"

# Record the machine + compiler + Minuit2 combo for reproducibility audit
{
    echo "# Reference data regen — machine fingerprint"
    echo "date_utc:    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "uname:       $(uname -a)"
    if [[ "$(uname)" == "Darwin" ]]; then
        echo "cpu:         $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
    fi
    echo "compiler:    $(${CXX:-c++} --version 2>/dev/null | head -1 || echo unknown)"
    echo "cmake:       $(cmake --version 2>/dev/null | head -1)"
    echo "minuit2_sha: $(cd "$REF_CPP" && git rev-parse HEAD 2>/dev/null)"
} > "$REF_DIR/_machine.txt"

echo ">>> Done. Reference data + machine fingerprint at $REF_DIR/"
ls -la "$REF_DIR/"

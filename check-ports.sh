#!/usr/bin/env bash
# Verify every buildable/runnable port against the C program's output.
#
# For each port whose toolchain is installed: build it if needed, run it
# at 100 and 1000 digits, and require byte-identical output files to the
# C program's. Ports without a toolchain are SKIPped. Exits nonzero if
# any port FAILs. Run via `make check` (which also verifies the C
# program itself against the reference digits first) or directly.
#
# Toolchain discovery covers PATH plus the user-local install locations
# this repository's tooling uses: ~/.dotnet, ~/.ghcup, ~/perl5 (GMP
# backend), opam, and ~/maple2024.

set -u
cd "$(dirname "$0")"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILURES=0

# Reference outputs from the C program.
make -s khinchin-fast
./khinchin-fast 100 "$WORK/ref100.txt" 2>/dev/null
./khinchin-fast 1000 "$WORK/ref1000.txt" 2>/dev/null

# check NAME -- CMD... ; CMD is run as "CMD digits outfile".
check() {
    local name="$1"; shift
    local ok=1
    for d in 100 1000; do
        if ! "$@" "$d" "$WORK/out.txt" >/dev/null 2>&1 \
            || ! cmp -s "$WORK/out.txt" "$WORK/ref$d.txt"; then
            ok=0
        fi
        rm -f "$WORK/out.txt"
    done
    if [ "$ok" = 1 ]; then
        printf '%-12s OK\n' "$name"
    else
        printf '%-12s FAIL\n' "$name"
        FAILURES=$((FAILURES + 1))
    fi
}

skip() { printf '%-12s SKIP (%s)\n' "$1" "$2"; }

# --- C++ ---
if command -v g++ >/dev/null; then
    g++ -O3 -march=native -std=c++17 -o "$WORK/khinchin-cpp" \
        ports/khinchin.cpp -lflint -lmpfr -lgmp -pthread 2>/dev/null \
        && check cpp "$WORK/khinchin-cpp" || skip cpp "build failed"
else skip cpp "no g++"; fi

# --- Rust ---
if command -v cargo >/dev/null; then
    cargo build --release --quiet \
        --manifest-path ports/khinchin-rs/Cargo.toml 2>/dev/null \
        && check rust ports/khinchin-rs/target/release/khinchin \
        || skip rust "build failed"
else skip rust "no cargo"; fi

# --- Julia ---
if command -v julia >/dev/null; then
    check julia julia -t auto ports/khinchin.jl
else skip julia "no julia"; fi

# --- Fortran ---
if command -v gfortran >/dev/null; then
    gfortran -O2 -fopenmp -o "$WORK/khinchin-f90" ports/khinchin.f90 \
        -lflint -lmpfr -lgmp 2>/dev/null \
        && check fortran "$WORK/khinchin-f90" || skip fortran "build failed"
else skip fortran "no gfortran"; fi

# --- PARI/GP ---
if command -v gp >/dev/null; then
    gp_run() {
        gp -q ports/khinchin.gp >/dev/null 2>&1 <<EOF
khinchin_write($1, "$2");
EOF
    }
    check gp gp_run
else skip gp "no gp"; fi

# --- Sage ---
if command -v sage >/dev/null; then
    check sage sage ports/khinchin_sage.sage
else skip sage "no sage"; fi

# --- Python (serial and two-region) ---
if command -v python3 >/dev/null && python3 -c 'import mpmath' 2>/dev/null; then
    check python python3 ports/khinchin.py
    check python-mt python3 ports/khinchin_mt.py
else skip python "no python3+mpmath"; fi

# --- Go ---
if command -v go >/dev/null; then
    check go go run ports/khinchin.go
else skip go "no go"; fi

# --- Java ---
if command -v java >/dev/null; then
    check java java ports/Khinchin.java
else skip java "no java"; fi

# --- C# ---
DOTNET=$(command -v dotnet || echo "$HOME/.dotnet/dotnet")
if [ -x "$DOTNET" ]; then
    export DOTNET_ROOT="${DOTNET_ROOT:-$(dirname "$DOTNET")}"
    "$DOTNET" build -c Release ports/khinchin-cs -v q >/dev/null 2>&1 \
        && check csharp ports/khinchin-cs/bin/Release/net9.0/khinchin-cs \
        || skip csharp "build failed"
else skip csharp "no dotnet"; fi

# --- Node.js ---
if command -v node >/dev/null; then
    check node node ports/khinchin.mjs
else skip node "no node"; fi

# --- Ruby ---
if command -v ruby >/dev/null; then
    check ruby ruby ports/khinchin.rb
else skip ruby "no ruby"; fi

# --- Perl (GMP backend picked up when present) ---
if command -v perl >/dev/null; then
    [ -d "$HOME/perl5/lib/perl5" ] \
        && export PERL5LIB="$HOME/perl5/lib/perl5${PERL5LIB:+:$PERL5LIB}"
    check perl perl ports/khinchin.pl
else skip perl "no perl"; fi

# --- Haskell ---
GHC=$(command -v ghc || echo "$HOME/.ghcup/bin/ghc")
if [ -x "$GHC" ]; then
    "$GHC" -O2 -o "$WORK/khinchin-hs" ports/khinchin.hs \
        -outputdir "$WORK/hs" >/dev/null 2>&1 \
        && check haskell "$WORK/khinchin-hs" || skip haskell "build failed"
else skip haskell "no ghc"; fi

# --- OCaml ---
if command -v opam >/dev/null || [ -x "$HOME/.local/bin/opam" ]; then
    PATH="$HOME/.local/bin:$PATH" opam exec -- ocamlfind ocamlopt \
        -package zarith -linkpkg ports/khinchin.ml \
        -o "$WORK/khinchin-ml" >/dev/null 2>&1 \
        && check ocaml "$WORK/khinchin-ml" || skip ocaml "build failed"
    rm -f ports/khinchin.cm* ports/khinchin.o
else skip ocaml "no opam"; fi

# --- Maple ---
MAPLE=$(command -v maple || echo "$HOME/maple2024/bin/maple")
if [ -x "$MAPLE" ]; then
    maple_run() {
        "$MAPLE" -q >/dev/null 2>&1 <<EOF
read "ports/khinchin.mpl":
khinchin_to_file($1, "$2"):
EOF
    }
    check maple maple_run
else skip maple "no maple"; fi

# --- Mathematica ---
if command -v wolframscript >/dev/null; then
    wl_run() {
        wolframscript -code \
            "Get[\"ports/khinchin.wl\"]; KhinchinToFile[$1, \"$2\"];" \
            >/dev/null 2>&1
    }
    check wolfram wl_run
else skip wolfram "no wolframscript"; fi

echo
if [ "$FAILURES" = 0 ]; then
    echo "All available ports byte-identical to the C program."
else
    echo "$FAILURES port(s) FAILED."
fi
exit "$FAILURES"

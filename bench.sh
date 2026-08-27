#!/usr/bin/env bash
# Regenerate the README's by-language speed table.
#
# Runs the C program and every port whose toolchain is installed at two
# digit counts (default 1000 and 10000) and prints a markdown table
# sorted by the larger measurement. All numbers are full wall-clock of
# one run each, C included — so interpreter/JIT startup is charged to
# every row equally (the README's curated table instead footnotes
# in-process numbers for Julia and Sage and uses best-of-N pairs; use
# repeated runs of this script when updating it). Output correctness is
# NOT checked here — run `make check` for that.
#
# Usage: ./bench.sh [digits_small digits_large]

set -u
cd "$(dirname "$0")"
D1=${1:-1000}
D2=${2:-10000}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
RESULTS="$WORK/results.tsv"
: > "$RESULTS"

wall() { # wall CMD... -> seconds
    local t0 t1
    t0=$(date +%s.%N)
    "$@" >/dev/null 2>&1
    t1=$(date +%s.%N)
    echo "$t0 $t1" | awk '{printf "%.2f", $2 - $1}'
}

row() { # row LABEL CMD... ; runs CMD D OUT for both digit counts
    local label="$1"; shift
    local s1 s2
    s1=$(wall "$@" "$D1" "$WORK/o.txt")
    s2=$(wall "$@" "$D2" "$WORK/o.txt")
    printf '%s\t%s\t%s\n' "$s2" "$label ${s1} s ${s2} s" "" >> "$RESULTS"
    printf '%-46s %8s s %10s s\n' "$label" "$s1" "$s2"
}

make -s khinchin-fast
row "C (khinchin-fast)" ./khinchin-fast
if [ -e "$HOME/opt/tuned-mathlibs/libgmp.so.10" ]; then
    tuned() { LD_PRELOAD="$HOME/opt/tuned-mathlibs/libgmp.so.10 $HOME/opt/tuned-mathlibs/libflint.so.22" ./khinchin-fast "$@"; }
    row "C + tuned libraries" tuned
fi

command -v g++ >/dev/null && g++ -O3 -march=native -std=c++17 \
    -o "$WORK/cpp" ports/khinchin.cpp -lflint -lmpfr -lgmp -pthread \
    2>/dev/null && row "C++ (ports/khinchin.cpp)" "$WORK/cpp"
command -v cargo >/dev/null && cargo build --release --quiet \
    --manifest-path ports/khinchin-rs/Cargo.toml 2>/dev/null \
    && row "Rust (ports/khinchin-rs)" ports/khinchin-rs/target/release/khinchin
command -v julia >/dev/null && row "Julia (ports/khinchin.jl)" julia -t auto ports/khinchin.jl
command -v gfortran >/dev/null && gfortran -O2 -fopenmp -o "$WORK/f90" \
    ports/khinchin.f90 -lflint -lmpfr -lgmp 2>/dev/null \
    && row "Fortran (ports/khinchin.f90)" "$WORK/f90"
if command -v gp >/dev/null; then
    gp_run() { gp -q ports/khinchin.gp >/dev/null 2>&1 <<EOF
khinchin_write($1, "$2");
EOF
    }
    row "PARI/GP (ports/khinchin.gp)" gp_run
fi
command -v sage >/dev/null && row "Sage (ports/khinchin_sage.sage)" sage ports/khinchin_sage.sage
if command -v python3 >/dev/null && python3 -c 'import mpmath' 2>/dev/null; then
    row "Python (ports/khinchin.py)" python3 ports/khinchin.py
    row "Python two-region (ports/khinchin_mt.py)" python3 ports/khinchin_mt.py
fi
if command -v go >/dev/null; then
    go build -o "$WORK/go-bin" ports/khinchin.go 2>/dev/null \
        && row "Go (ports/khinchin.go)" "$WORK/go-bin"
fi
command -v java >/dev/null && row "Java (ports/Khinchin.java)" java ports/Khinchin.java
DOTNET=$(command -v dotnet || echo "$HOME/.dotnet/dotnet")
if [ -x "$DOTNET" ]; then
    export DOTNET_ROOT="${DOTNET_ROOT:-$(dirname "$DOTNET")}"
    "$DOTNET" build -c Release ports/khinchin-cs -v q >/dev/null 2>&1 \
        && row "C# (ports/khinchin-cs)" ports/khinchin-cs/bin/Release/net9.0/khinchin-cs
fi
command -v node >/dev/null && row "Node.js (ports/khinchin.mjs)" node ports/khinchin.mjs
command -v ruby >/dev/null && row "Ruby (ports/khinchin.rb)" ruby ports/khinchin.rb
if command -v perl >/dev/null; then
    [ -d "$HOME/perl5/lib/perl5" ] \
        && export PERL5LIB="$HOME/perl5/lib/perl5${PERL5LIB:+:$PERL5LIB}"
    row "Perl (ports/khinchin.pl)" perl ports/khinchin.pl
fi
GHC=$(command -v ghc || echo "$HOME/.ghcup/bin/ghc")
[ -x "$GHC" ] && "$GHC" -O2 -o "$WORK/hs" ports/khinchin.hs \
    -outputdir "$WORK/hsbuild" >/dev/null 2>&1 \
    && row "Haskell (ports/khinchin.hs)" "$WORK/hs"
if command -v opam >/dev/null || [ -x "$HOME/.local/bin/opam" ]; then
    PATH="$HOME/.local/bin:$PATH" opam exec -- ocamlfind ocamlopt \
        -package zarith -linkpkg ports/khinchin.ml -o "$WORK/ml" \
        >/dev/null 2>&1 && row "OCaml (ports/khinchin.ml)" "$WORK/ml"
    rm -f ports/khinchin.cm* ports/khinchin.o
fi
MAPLE=$(command -v maple || echo "$HOME/maple2024/bin/maple")
if [ -x "$MAPLE" ]; then
    maple_run() { "$MAPLE" -q >/dev/null 2>&1 <<EOF
read "ports/khinchin.mpl":
khinchin_to_file($1, "$2"):
EOF
    }
    row "Maple (ports/khinchin.mpl)" maple_run
fi
if command -v wolframscript >/dev/null; then
    wl_run() { wolframscript -code \
        "Get[\"ports/khinchin.wl\"]; KhinchinToFile[$1, \"$2\"];" \
        >/dev/null 2>&1; }
    row "Mathematica (ports/khinchin.wl)" wl_run
fi

echo
echo "Markdown, sorted by the ${D2}-digit column:"
echo
echo "| Implementation | ${D1} digits | ${D2} digits |"
echo "|---|---:|---:|"
sort -g "$RESULTS" | while IFS=$'\t' read -r _ line _; do
    # line = "label s1 s s2 s"
    set -- $line
    label=""
    while [ $# -gt 4 ]; do label="$label $1"; shift; done
    echo "|${label} | $1 $2 | $3 $4 |"
done

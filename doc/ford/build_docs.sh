#!/usr/bin/env bash
# Generates the FORD API documentation for AGILE.
#
# FORD parses plain Fortran, not fypp (`#:if`, `${...}$`, ...), so this
# script first expands every `.fpp` file under src/ with fypp -- the same
# tool and flags the build system uses (see make/30-fypp.mk) -- into
# `.f90`, then points FORD at the result.
#
# AGILE's physics modules (hd/mhd/ffhd/srhd) are mutually-exclusive fypp
# (`#:if PHYS == '...'`) branches, so a single fypp pass with one PHYS
# value only shows that physics module's variant of any file whose content
# depends on it. To document all four physics modules in one site:
#  1. every `.fpp` file is preprocessed once per physics module (using
#     docgen_{mhd,hd,ffhd,srhd}.par), into four staging directories.
#  2. merge_physics_variants.py compares each file's four variants: files
#     byte-identical across all four (the vast majority -- anything with
#     no PHYS-dependent content) are kept once under their real name;
#     files that differ are kept once per physics module, suffixed
#     `_<phys>` (filename and any `module`/`end module` statement), so all
#     four variants can coexist in one FORD site without name collisions.
#     This auto-discovers which files are physics-sensitive -- no
#     hand-maintained file list.
#
# Suffixed module names (e.g. `mod_physics_hd`) are a documentation-only
# convenience -- they aren't real compiled Fortran identifiers; the real,
# single module (e.g. `mod_physics`) compiles to whichever PHYS was
# selected at actual build time.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
agile_dir="$(cd "$script_dir/../.." && pwd)"
stage_dir="$script_dir/.stage"
src_out="$script_dir/src"
html_out="$script_dir/html"

if ! command -v fypp >/dev/null 2>&1; then
    echo "error: fypp not found. Run 'uv sync --extra docs' and activate .venv first." >&2
    exit 1
fi
if ! command -v ford >/dev/null 2>&1; then
    echo "error: ford not found. Run 'uv sync --extra docs' and activate .venv first." >&2
    exit 1
fi

# Prints fypp flags (one per line) derived from a docgen_*.par file.
fypp_flags_for() {
    python "$agile_dir/make/config_reader.py" < "$1" \
        | grep '^fypp_flags += ' \
        | sed -e 's/^fypp_flags += //' -e "s/\\\\'/'/g"
}

rm -rf "$stage_dir" "$src_out" "$html_out"

for phys in mhd hd ffhd srhd; do
    echo "=== Preprocessing pass (PHYS=$phys) ==="
    out_dir="$stage_dir/$phys"
    mkdir -p "$out_dir"

    fypp_flags=(-M "$agile_dir/make")
    while read -r flag; do
        fypp_flags+=("$flag")
    done < <(fypp_flags_for "$script_dir/docgen_$phys.par")

    while IFS= read -r -d '' f; do
        out="$out_dir/$(basename "${f%.fpp}").f90"
        fypp "${fypp_flags[@]}" "$f" "$out"
    done < <(find "$agile_dir/src" -name '*.fpp' -print0)
done

echo "Merging physics variants -> $src_out"
mkdir -p "$src_out"
python "$script_dir/merge_physics_variants.py" "$stage_dir" "$src_out"
rm -rf "$stage_dir"

echo "Rewriting same-line trailing doc comments (see fix_inline_docs.py)"
python "$script_dir/fix_inline_docs.py" "$src_out"

echo "Running FORD"
cd "$script_dir"
ford project.md

echo "Documentation generated in $html_out/index.html"

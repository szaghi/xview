# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

xview is a single-executable CLI that post-processes the output of the Xall/Xnavis CFD solver (overset/Chimera
structured multi-block grids). It reads grid (`.grd`), topology (`.icc`) and restart (`rst`) files, or
extracted-subzone (`esz`) files, computes derived fields, and writes Tecplot (`.plt` via the bundled `tecio`
library) and/or VTK (`.vts` per block plus a `.vtm` multiblock index) files.

## Build

Build with FoBiS; the `fobos` file lives at the repo root:

```bash
fobis fetch                          # clone the dependencies into src/third_party/ (once; --update to refresh)
fobis build --mode xview-gnu-debug   # bin/xview; also: xview-gnu, xview-intel, xview-intel-debug, xview-nvf
fobis build --lmodes                 # list modes
fobis rule --ex clean                # also deletes bin/, exe/ and the built docs
fobis rule --ex makedoc              # formal API + VitePress site in docs/.vitepress/dist
```

- **Dependencies:** they come from `fobos [dependencies]` via `fobis fetch`, with `use=sources`. Each clone is compiled
  inline, so its own nested `src/third_party` stays empty. The list is therefore the **flat transitive closure**:
  PENF, FACE, BeFoR64, StringiFor, FLAP, VecFor, FoXy, VTKFortran.
  - `$EXDIRS` excludes each clone's `src/tests/ src/third_party/ scripts/ docs/`, and also `FoXy/benchmarks/`.
  - When adding a dependency, add its own dependencies and its `$EXDIRS` entries too. Then check that
    `find src/third_party -name '*.[fF]90'` shows nothing outside `*/src/lib/`.
  - `src/third_party/` is git-ignored, except its `.gitignore` and `.deps_config.ini`. `fobos.lock` lands there and is
    not tracked.
- **Output:** the application builds into `bin/xview`, and the tests build into `exe/` (the scaffold layout).
  - The vendored `lib/tecio/libtecio.so` is linked with `-Wl,-rpath,'$ORIGIN/../lib/tecio'` (`$TECIO_RPATH`), so
    binaries in `bin/` and `exe/` find it without `LD_LIBRARY_PATH`.
  - Do **not** add `stdc++` to `$TECIO_LIBS`. With `lib_dir = lib/tecio/`, an old vendored `libstdc++` there once made
    `-lstdc++` bind the executable to the obsolete `libstdc++.so.5`. `libtecio.so` itself needs only the system
    `libstdc++.so.6`.
- **nvfortran and VecFor operators:** nvfortran (every version from 25.11 to 26.9) mis-compiles VecFor's defined
  operators `.dot.`, `.cross.` and `.ortho.`.
  - Symptoms: an internal compiler error ("mkexpr1: bad const"), or a misparse such as reading
    `a .dot. self%NFi(0,j,k)` as `(a .dot. self)%NFi(...)`. Which operand forms trigger it is erratic.
  - The code therefore uses the type-bound calls: `a%dotproduct(b)`, `a%crossproduct(b)`, `a%orthogonal(b)`.
    **Never write these operators in xview sources.**
  - When the left operand is an expression, store it in a local first (as `dm` in `compute_yplus`).
  - Keep the operator's grouping. A defined binary operator has the lowest precedence, so `a.dot.b/c` means
    `a%dotproduct(b/c)`, not `a%dotproduct(b)/c`.
  - Arithmetic operators (`+ - * /`) on vectors compile fine.
- **Stale README help:** the CLI help in the README uses the old single-dash switches. Use `xview --help`.

## Running / tests

```bash
fobis build --mode tests-gnu-debug && bash scripts/run_tests.sh   # what CI runs, without coverage
fobis rule --ex makecoverage-analysis                              # CI entry point: tests + coverage + reports
```

- **How tests work:** each `src/tests/xview_test_*.F90` is a self-checking program. It loads the
  `src/tests/xnavis-level-set/input/` case through the library API and fails with `error stop`.
  `scripts/run_tests.sh` (scaffold-managed, do not edit) runs every executable in `exe/` with no arguments, **from
  the repository root**, so tests open their input data by repository-relative paths.
- **Mirror the CLI:** when a test needs gradient fields (vorticity, lambda2, liutex, grad_p, div2LT…), pass
  `compute_metrics=.true.` to `file_grd%load_file` and `correct_metrics=.true.` to `file_icc%load_file`, as `parse_cli`
  does implicitly. Without them the metric arrays are not allocated, and the gradient walks out of bounds.
- **Smoke scripts:** the shell scripts under `src/tests/*/` are runs of `bin/xview` that check nothing.
  For example: `cd src/tests/xnavis-level-set && ./xview_postprocess_xnavis_level_set_test.sh`.

That script also shows the canonical invocations (single file, grd+icc, grd+icc+rst tuple, `--glob`, `--patch`).
`utils/run_xview.sh` is the production driver script that wraps the CLI.

CLI shape: global switches go **before** the command, and command-specific switches go after it:

```
xview [--grd ..] [--icc ..] [--rst ..] [--ipath ..] [--opath ..] [--glob] [--procinput|--blocks-map ..] postprocess [--cell] [--tec yes|no] [--vtk yes|no] [--lambda2] ...
```

FLAP defines the commands `postprocess`, `interpolate` and `aeroacustic`. The program body (`src/app/xview.F90`)
dispatches only `postprocess`. The other two have CLI switches but no implementation.

## Architecture

All library code is in `src/lib/`, one module per file. The naming convention is:

- `xview_file_*_object`: a file-level reader/writer. It extends `file_object` (which holds `filename`,
  `is_loaded` and `is_file_present`) and holds an array of blocks.
- `xview_block_*_object`: per-block data. `block_grd` holds nodes and metrics, `block_icc` holds the
  cell-connectivity/BC flags, `block_rst` holds the solution plus **all derived fields**, and `block_esz`
  aggregates the grd, icc and rst data of one subzone.
- `xview_file_tec_object` / `xview_file_vtk_object` are the output writers. `file_ptc` holds patch extents
  (`patches.dat`) and writes the aeroacoustic patch output. `file_procinput`, `file_blksmap` and `file_mbpar` are
  auxiliary inputs that map files and ranks to global block numbers.
- `xview_ui_object` wraps FLAP: `set_cli` declares the switches and `parse_cli` copies them into flat `ui%...`
  fields (`compute_*`, `is_*`, paths).

Data flow in `src/app/xview.F90`:
1. `get_filenames` resolves the single files, or with `--glob` the per-process file lists (`*.grd.p???`, rst files
   by `--sname-rst` and timestep for esz).
2. `postprocess_gis` (grd/icc/rst) and `postprocess_esz` (extracted subzones, looping over timesteps × blocks) load
   the files, then call the VTK and/or Tecplot writers for each block.
3. The `.vtm` file that groups all `.vts` files is written at the end. `files_blocks(3, nfiles)` threads the
   grouping through.

Derived fields (lambda2, qfactor, liutex, helicity, vorticity, grad_p, Lighthill tensor LT/div2LT, k-ratio,
yplus, tau, div_tau, loads) are computed in `xview_block_rst_object.F90`, mainly `compute_aux` and its helpers.
Each field is an `allocatable` component gated by a `has_<field>` flag. Allocation happens in `alloc`, with ghost
cells `gc(1:6)`.

**Adding a derived field touches all of these places, and missing one fails silently:**
1. `ui_object`: a switch in `set_cli`, the flag in `parse_cli`, and a line in `description`.
2. The `compute_<field>` optional arguments, threaded through `file_rst%load_file` → `block_rst%load_solution` →
   `init`, and through `file_esz%load_file` → `block_esz%load_solution`.
3. Both `load_file` call sites in `xview.F90`.
4. `block_rst`: the component, the `has_` flag, `alloc`, `destroy`, the computation, and `interpolate_at_nodes`.
5. The writers: `file_vtk` and `file_tec` both write a field when its array is `allocated`. In `file_tec` the field
   needs two edits: a `tec_var` call in `save_block_file_tec` and an `add_var` call in `initialize`, at the **same
   position**. The header names come from that list. If the two orders differ, the Tecplot columns are silently
   shifted. To verify, check that each zone in an `--ascii` output holds as many values as the header declares
   variables.

## Conventions specific to this repo

- **Scaffold:** the project kit (CI and docs workflows, `scripts/`, `cliff.toml`, `docs/package.json`, the
  contributing page) is managed by `fobis scaffold`. Check it with `fobis scaffold status` and update it with
  `fobis scaffold sync`.
  - Edit those files upstream, in the scaffold templates, not here.
  - `fpm.toml` is skipped on purpose (`fobos [scaffold] skip`): fpm cannot link the vendored tecio.
  - `.github/workflows/ci.yml` and `docs.yml` are disabled for now, by being listed in `[scaffold] skip`. Nothing runs
    on a branch push. `release.yml` still runs on a `vX.Y.Z` tag push.
  - `docs/ford.md`, `docs/.vitepress/config.mts` and the `docs/guide/` pages are project-owned.
- **Releases:** `scripts/release.sh vX.Y.Z` regenerates `CHANGELOG.md` with git-cliff from the Conventional Commits.

- `.F90` files use the `-cpp` preprocessor and `-frealloc-lhs` (GNU). Kinds come from PENF (`I4P`, `R8P`, ...),
  strings are StringiFor `string`, and 3-vectors are VecFor `vector`.
- Grid and restart files are Fortran `unformatted` sequential, written by Xnavis, in native endianness. No
  `convert=` is set.
- `building-errors.log` at the root is a scratch log, not a build artifact the code depends on.

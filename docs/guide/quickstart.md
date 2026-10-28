---
title: Quick Start
---

# Quick Start

Global switches (input files and paths) go **before** the command, the command switches after it:

```bash
xview [global switches] postprocess [postprocess switches]
```

`xview --help` lists all the switches; `xview postprocess --help` lists the ones of the command.

## One grid/topology/restart tuple

```bash
xview --blocks-map blks_map.ini --myrank 2 \
      --grd CC/cc.02.grd.p002 --icc CC/cc.02.p002 --rst RST/boa_00.02.p002 \
      --ipath input/ --opath output/ \
      postprocess --cell
```

`--cell` saves the variables at cell centers instead of interpolating them at the nodes. `--tec yes|no` and
`--vtk yes|no` select the output formats, `--ascii` writes ascii instead of binary files.

## All the files of a parallel run

With `--glob`, the file names become base names and xview processes every `*.p???` file found:

```bash
xview --blocks-map blks_map.ini --glob --grd CC/cc --icc CC/cc --rst RST/boa \
      --ipath input/ --opath output/ \
      postprocess --cell
```

## Derived fields

Add the fields to the `postprocess` command, for example:

```bash
xview ... postprocess --cell --lambda2 --qfactor --liutex --vorticity --LT --div2LT
```

Vortex identification: `--lambda2`, `--qfactor`, `--liutex`, `--helicity`, `--vorticity`. Aeroacoustics:
`--LT` (Lighthill tensor), `--div2LT` (its double divergence). Others: `--gradp`, `--k-ratio`, and on wall patches
(`--patch`) `--yplus`, `--tau`, `--div-tau`, `--loads`.

The test case `src/tests/xnavis-level-set/` holds a small level-set run with the command lines of the examples above.

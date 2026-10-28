---
title: Installation
---

# Installation

## Requirements

- A Fortran 2018 compiler: GNU gfortran (tested with 14.3 and 16.2), NVIDIA nvfortran (tested with 26.9) or Intel
  `ifx` (build modes present, not yet verified).
- [FoBiS](https://github.com/szaghi/FoBiS) 3.8+ to fetch the dependencies and build.
- Linux x86-64: the bundled Tecplot I/O library `lib/tecio/libtecio.so` is a Linux x86-64 binary.

## Build

```bash
git clone https://github.com/szaghi/xview
cd xview
fobis fetch                     # clone the dependencies into src/third_party/
fobis build --mode xview-gnu    # build bin/xview
```

`fobis build --lmodes` lists the other modes (`xview-gnu-debug`, `xview-intel`, `xview-intel-debug`, `xview-nvf`).
The executable finds `libtecio.so` through its run path, so it needs no `LD_LIBRARY_PATH` as long as it stays in
`bin/` next to `lib/`.

## Dependencies

`fobis fetch` clones them from the `[dependencies]` section of `fobos`:
[PENF](https://github.com/szaghi/PENF), [FACE](https://github.com/szaghi/FACE),
[BeFoR64](https://github.com/szaghi/BeFoR64), [StringiFor](https://github.com/szaghi/StringiFor),
[FLAP](https://github.com/szaghi/FLAP), [VecFor](https://github.com/szaghi/VecFor),
[FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy) and [VTKFortran](https://github.com/szaghi/VTKFortran).
`fobis fetch --update` refreshes them.

## Tests

```bash
fobis build --mode tests-gnu
bash scripts/run_tests.sh
```

---
title: About
---

# About

xview is a command line tool that post-processes the output of the Xall/Xnavis CFD solver, which works on overset
(Chimera) structured multi-block grids. It reads:

- grid files (`.grd`), with or without ghost cells;
- topology files (`.icc`), the cell-connectivity and boundary-condition flags;
- restart files, the solution;
- extracted-subzone files of unsteady runs.

It computes derived fields on demand and writes Tecplot files (`.plt` binary through the bundled `tecio` library, or
ascii `.dat`) and VTK files (one `.vts` per block plus a `.vtm` multiblock index).

## Authors

- Stefano Zaghi
- Cristiano Andolfi
- Andrea Di Mascio

## License

xview is distributed under the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.html).

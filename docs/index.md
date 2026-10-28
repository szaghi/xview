---
layout: home

hero:
  name: xview
  text: post-processing for Xall/Xnavis
  tagline: Load Xall/Xnavis grid, topology and restart files, compute derived fields, and write Tecplot or VTK files for ParaView and Tecplot.
  actions:
    - theme: brand
      text: Quick Start
      link: /guide/quickstart
    - theme: alt
      text: Installation
      link: /guide/installation
    - theme: alt
      text: API
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/xview

features:
  - title: Xall/Xnavis inputs
    details: Grid (.grd), topology (.icc) and restart files of overset multi-block grids, single tuples or whole parallel runs found by glob; extracted subzones of unsteady runs.
  - title: Derived fields
    details: Vortex identification (lambda2, Q-criterion, Liutex, helicity, vorticity), pressure gradient, Lighthill tensor and its double divergence, y+, wall shear stress, loads on patches.
  - title: Tecplot and VTK
    details: Binary or ascii Tecplot files through the tecio library; VTK structured blocks grouped in a multiblock .vtm file for ParaView.
---

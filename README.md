<a name="top"></a>

# XVIEW

> xview, post-processing and analysis tool for Xall/Xnavis CFD solver

### Authors

+ Stefano Zaghi, [stefano.zaghi@cnr.it](stefano.zaghi@cnr.it)
+ Cristiano Andolfi, [cristiando93@gmail.com](cristiando93@gmail.com)
+ Andrea di Mascio, [andrea.dimascio@univaq.it](andrea.dimascio@univaq.it)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()

[![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]()
[![Build Status](https://travis-ci.org/szaghi/xview.svg?branch=master)](https://travis-ci.org/szaghi/xview)
[![Coverage Status](https://img.shields.io/codecov/c/github/szaghi/xview.svg)](http://codecov.io/github/szaghi/xview?branch=master)

#### <a name="toc">Table of Contents

* [What is xview?](#what)
* [Todos](#todos)
* [Requirements](#requirements)
* [Copyrights](#copyrights)
* [Download xview](#download)
* [Compiling Instructions](#compile)
* [Usage](#usage)
  + [Main Help](#usage-help)
  + [Post-processing only mesh files](#usage-only-mesh)
  + [Post-processing mesh and solutions files](#usage-mesh-sol)
  + [Utilities](#utilities)

Go to [Top](#top) or [Toc](#toc)

#### Issues

[![GitHub issues](https://img.shields.io/github/issues/szaghi/xview.svg)]()
[![Ready in backlog](https://badge.waffle.io/szaghi/xview.png?label=ready&title=Ready)](https://waffle.io/szaghi/xview)
[![In Progress](https://badge.waffle.io/szaghi/xview.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/xview)
[![Open bugs](https://badge.waffle.io/szaghi/xview.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/xview)

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-v9.2.0+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/Intel-v2024+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/g95-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/NAG-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/PGI-not%20tested-yellow.svg)]()

## <a name="what"></a>What is xview?

xview loads Xall/Xnavis outputs and produces post-processed files.

Go to [Top](#top) or [Toc](#toc)

## <a name="todos"></a>Todos

+ any feature request is welcome!

Go to [Top](#top) or [Toc](#toc)

## <a name="download"></a>Download xview

If you use `git` it could be convenient to clone this repository:
```bash
git clone --recursive https://github.com/szaghi/xview
```
Other 2 possibilities are:

1. use the GitHub **Download ZIP** button on the right sidebar of this page;
2. download one of the releases in the [release page](https://github.com/szaghi/xview/releases), also listed at the end of this page.

Go to [Top](#top) or [Toc](#toc)

## <a name="requirements"></a>Requirements

+ Modern Fortran Compiler (standard 2003+);
+ a lot of patience with the authors.

xview is developed on a GNU/Linux architecture. For Windows architecture there is no support, however it should be work out-of-the-box.

Go to [Top](#top) or [Toc](#toc)

## <a name="Copyrights"></a>Copyrights

xview is an open source project, it is distributed under the [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest to use, to develop or to contribute to xview is welcome.

Go to [Top](#top) or [Toc](#toc)

## <a name="compile"></a>Compiling Instructions

xview has been developed on GNU/Linux architectures. Other OS are not supported (and in general there is no best alternative to GNU/Linux :-).

xview have been successfully compiled with the following compilers:

+ GNU gfortran (version 4.7.0 or higher);
+ Intel Fortran Compiler ifort (version 12.0 or higher)

xview is constituted by several modules. Therefore there are many dependences. The most easy way to compile the code is to start with the provided makefile thus it is necessary that the system has "Make" program (preferably GNU make http://www.gnu.org/software/make/).

The provided makefile has several options. It has one rule that prints all options available and the default settings. Typing in the shell prompt: `code make help` the following output will be printed:

```shell
 Make options of xview code

 Compiler choice: COMPILER=intel => default
  COMPILER=gnu   => GNU gfortran
  COMPILER=intel => Intel Fortran

 Compiling options
  DEBUG=yes(no)    => on(off) debug                  (default no)
  F03STD=yes(no)   => on(off) check standard fortran (default no)
  OPTIMIZE=yes(no) => on(off) optimization           (default no)
  OPENMP=yes(no)   => on(off) OpenMP directives      (default no)
  BIGEIN=yes(no)   => on(off) Big Endian input files (default no)

 Preprocessing options
  R16P=yes(no) => on(off) definition of real with "128 bit" (default no)

 Executable directory
  DEXE="your_path" => directory where exe is placed (default ~/bin/)

 External libraries
  TECIO=yes(no) => on(off) Tecplot IO library linking (default )

 Provided Rules
  Defualt rule => ~/bin/xview
  help         => printing this help message
  cleanobj     => cleaning compiled object
  cleanmod     => cleaning .mod files
  cleanmsg     => cleaning make-log massage files
  cleanexe     => cleaning executable files
  clean        => running cleanobj, cleanmod and cleanmsg
  cleanall     => running clean and cleanexe
  tar          => creating a tar archive of the project
```

For example compiling in debug mode with the Intel Fortran compiler you can type:

```shell
make DEBUG=yes COMPILER=intel
```

Go to [Top](#top) or [Toc](#toc)

## <a name="usage"></a>Usage

### <a name="usage-help">Main Help

xview is is a Command Line Tool. To list the available options run it as following:

```shell
./xview -h
```

this help will echoed

```shell
xview, post-processor for Xnavis
usage: xview  -g value [-o value] [-i value] [-s value] [-ngc] [-cell] [-ls] [-nt] [-eq value] [-vordet] [-ascii] [-tec value] [-vtk value] [-proc value] [-os value] [-vb] [--help] [--version]
Required switches:
   -g value
          Grid file (.grd)
Optional switches:
   -o value
          default value unset
          output file name; default is basename of grd file with the proper extension
   -i value
          default value unset
          ICC file
   -s value
          default value unset
          solution file name; if passed the solution variables are saved
   -ngc
          default value .false.
          mesh without ghosts cells, as geogrd output
   -cell
          default value .false.
          variables other than nodes coord. are cell centered
   -ls
          default value .false.
          solution with level set
   -nt
          default value .false.
          no turbulent model used
   -eq value, value in: (0,1,2)
          default value 1
          equations turbulent model
   -vordet
          default value .false.
          computing variables for vortices identification
   -ascii
          default value .false.
          write ascii output files
   -tec value, value in: (yes,no)
          default value yes
          write output Tecplot files
   -vtk value, value in: (yes,no)
          default value no
          write output VTK files
   -proc value
          default value -1
          process number for global block numeration if proc.input is found
   -os value, value in: (UIX,WIN)
          default value UIX
          type of Operating System
   -vb
          default value .false.
          Verbose output
   --help, -h
          Print this help message
   --version, -v
          Print version

Examples:
   xview -g xship.grd -o grid
   xview -g cc.01.grd -i cc.01 -o mesh.01
   xview -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01
```

### <a name="usage-only-mesh">Post-processing only mesh files

#### GRD files (no ghost cells)

```shell
  xview -g xship.grd -o mesh
```

A file named `mesh.plt` is generated.

#### Overset/Overott files (with ghost cells)

```shell
  xview -g cc.01.grd -i cc.01 -o mesh.01
```

A file named `mesh.01.plt` is generated.

### <a name="usage-mesh-sol">Post-processing mesh and solutions files

```shell
  xview -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01
```

A file named `sol.01.plt` is generated.

### <a name="utilities">Utilities

To be written.

Go to [Top](#top) or [Toc](#toc)

### Features

+ KISS, keep it simple and stupid;
+ simulate different types of flows;
+ use high order finite difference scheme;
+ exploit Adaptive Mesh Refinement (AMR) to accurate simulate complex geometries;
+ exploit Immersed Boundary (IB) techniques to easy handle complex moving geometries;
+ ready for very High Performance Computing (HPC) by means of GPU parallel paradigm exploitation;

## Contributing

The project home is hosted on GitHub and the development is organized by means of multiple branches as described in the following.

### Branches Organization

Currently, ADAM development tree contains the following branches:

+ master
+ develop
+ fortran-cuda-stable
+ fortran-cuda-develop
+ openmp-gpu-develop
+ euler-cpu-develop

#### master

Soon or later a stable v1.0.0 will born and `master` branch will be its home, for now it is an empty branch.

#### develop

It has been the home for the Fortran-CUDA development until now. In 15th September 2023 it has been copied into branch `fortran-cuda-stable` and it is currently freeze.

#### fortran-cuda-stable

It is a **reference** branch, used for comparison reason, a sort of *named commit*.

#### fortran-cuda-develop

It contains the ongoing Fortran-CUDA development.

Currently, the main CUDA development is focused into the NASTO (CUDA) application. More details about NASTO app can be found in its [readme](src/app/nasto/README.md).

#### openmp-gpu-develop

It contains the ongoing OpenMP development, GPU offloading.

#### euler-cpu-develop

It contains the ongoing development of a very simplified application for Euler flows simulation with a simple CPU backend.

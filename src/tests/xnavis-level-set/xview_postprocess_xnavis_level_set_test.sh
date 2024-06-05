#!/bin/bash

# test single grd postprocessing
./xview --blocks-map blks_map.ini --myrank 1 --grd CC/cc.02.grd.p001 --ipath input/ --opath output-g/ --verbose postprocess --cell
# test single grd/icc postprocessing
./xview --blocks-map blks_map.ini --myrank 3 --grd CC/cc.02.grd.p003 --icc CC/cc.02.p003 --ipath input/ --opath output-gi/ --verbose postprocess --cell
# test single tuple (Grd/Icc/rSt, gis) of files postprocessing
./xview --blocks-map blks_map.ini --myrank 2 --grd CC/cc.02.grd.p002 --icc CC/cc.02.p002 --rst RST/boa_00.02.p002 --ipath input/ --opath output-gis/ --verbose postprocess --cell
# test glob files search postprocessing
./xview --blocks-map blks_map.ini --glob --grd CC/cc --icc CC/cc --rst RST/boa --ipath input/ --opath output-glob/ --verbose postprocess --cell

# test single tuple (Grd/Icc/rSt, gis) of files postprocessing with patches
./xview --blocks-map blks_map.ini --myrank 2 --grd CC/cc.02.grd.p002 --icc CC/cc.02.p002 --rst RST/boa_00.02.p002 --ipath input/ --opath output-gis-patch-01/ --verbose postprocess --cell --patch 1
# test glob files search postprocessing
./xview --blocks-map blks_map.ini --glob --grd CC/cc --icc CC/cc --rst RST/boa --ipath input/ --opath output-glob-patch-01/ --verbose postprocess --cell --patch 1

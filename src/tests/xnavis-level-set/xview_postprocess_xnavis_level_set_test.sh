#!/bin/bash

# test glob files search postprocessing
./xview --blocks-map blks_map.ini --glob --grd CC/cc --icc CC/cc --rst RST/boa --ipath input/ --verbose postprocess
# test single file search postprocessing
./xview --blocks-map blks_map.ini --myrank 2 --grd CC/cc.02.grd.p002 --icc CC/cc.02.p002 --rst RST/boa_00.02.p002 --ipath input/ --verbose postprocess

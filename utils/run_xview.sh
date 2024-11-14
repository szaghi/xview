#!/bin/bash

# options

# verbose output
VERBOSE='no'
# grid level
LEVEL=1
# grd basename
GRD='CC/cc'
# icc basename
ICC='CC/cc'
# solution basename
SOL='RST/kde'
# output path
OPATH='postprocessing/'
# output basename
ONAME='kde'
# postprocess volume
VOLUME='yes'
# postprocess wall surfaces
WALL='no'
# cell-centered
CELL='--cell'
# compute auxiliary varibles in postprocess mode, default no
LAMBDA2=''
QFACTOR=''
HELICITY=''
VORTICITY=''
GRADP=''
DIV2LT=''
KRATIO=''
YPLUS=''
TAU=''
DIVTAU=''

#parsing command line
while [ $# -gt 0 ]; do
  case "$1" in
    "-verbose")
      shift; VERBOSE='yes'
      ;;
    "-l")
      shift; LEVEL=$1
      ;;
    "-grd")
      shift; GRD=$1
      ;;
    "-icc")
      shift; ICC=$1
      ;;
    "-sol")
      shift; SOL=$1
      ;;
    "-opath")
      shift; OPATH=$1
      ;;
    "-oname")
      shift; ONAME=$1
      ;;
    "-volume")
      shift; VOLUME=$1
      ;;
    "-wall")
      shift; WALL=$1
      ;;
    "-node")
      shift; CELL=''
      ;;
    "-lambda2")
      shift; LAMBDA2='--lambda2'
      ;;
    "-qfactor")
      shift; QFACTOR='--qfactor'
      ;;
    "-helicity")
      shift; HELICITY='--helicity'
      ;;
    "-vorticity")
      shift; VORTICITY='--vorticity'
      ;;
    "-gradp")
      shift; GRADP='--gradp'
      ;;
    "-div2LT")
      shift; DIV2LT='--div2LT'
      ;;
    "-k-ratio")
      shift; KRATIO='--k-ratio'
      ;;
    "-yplus")
      shift; YPLUS='--yplus'
      ;;
    "-tau")
      shift; TAU='--tau'
      ;;
    "-div-tau")
      shift; DIVTAU='--div-tau'
      ;;
    *)
      echo; echo "Unknown switch $1"; exit 1
      ;;
  esac
  shift
done

RE=$(grep -i Reyn mb.par | awk '{print $1}')

if [ "$VOLUME" == "yes" ] ; then
   xview --blocks-map blk_groups.ini --glob --grid-level $LEVEL --grd $GRD --icc $ICC --rst $SOL --opath $OPATH --out $ONAME --verbose postprocess $CELL $LAMBDA2 $QFACTOR $HELICITY $VORTICITY $GRADP $DIV2LT $KRATIO $YPLUS $TAU $DIVTAU
   if [ "$VERBOSE" == "yes" ] ; then
     echo "xview --blocks-map blk_groups.ini --glob --grid-level $LEVEL --grd $GRD --icc $ICC --rst $SOL --opath $OPATH --out $ONAME --verbose postprocess $CELL $LAMBDA2 $QFACTOR $HELICITY $VORTICITY $GRADP $DIV2LT $KRATIO $YPLUS $TAU $DIVTAU"
   fi
fi
# if [ "$WALL" == "yes" ] ; then
# fi

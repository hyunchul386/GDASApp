#!/bin/bash
set -ex

srcdir=$1

here=`pwd`
echo "here : " $here
echo "data : " $data

if [ ! -f "${data}/ocn_da_2021_03_22_08.nc ufs_input.nc" ]; then
#---- only test
   cp /scratch1/NCEPDEV/stmp4/Hyun-Chul.Lee/GW-gdasapp_test/ocn_da_2021_03_22_08.nc ${data}/.
fi
if [ ! -e "ufs_input.nc" ]; then
   ln -s ${data}/ocn_da_2021_03_22_08.nc ufs_input.nc
fi
#--- input yaml
if [ ! -f "./testinput/plot_front.yaml" ]; then
  cp ${srcdir}/test/soca/testinput/plot_front.yaml ./testinput/.
fi
python ${srcdir}/ush/soca/plot_front_ufs.py
re=$?
export err=$re

if [ ! -d testoutput ]; then
   mkdir testoutput
fi
mv front_output*png testoutput/.
exit $err


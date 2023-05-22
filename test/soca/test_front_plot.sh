#!/bin/bash
set -ex

srcdir=$1

here=`pwd`
echo "here : " $here
echo "data : " $data
#module purge
#module use  ${srcdir}/sorc/gdas.cd/modulefiles
#module load EVA/hera

if [ ! -f "${data}/ocn_da_2021_03_22_08.nc ufs_input.nc" ]; then
#---- onlt test
   cp /scratch1/NCEPDEV/stmp4/Hyun-Chul.Lee/GW-gdasapp_test/ocn_da_2021_03_22_08.nc ${data}/.
fi
if [ ! -e "ufs_input.nc" ]; then
   ln -s ${data}/ocn_da_2021_03_22_08.nc ufs_input.nc
fi
python ${srcdir}/ush/soca/plot_front_ufs.py
re=$?
export err=$re

if [ ! -d testoutput ]; then
   mkdir testoutput
fi
mv front_output*png testoutput/.
exit $err

#module purge
#module use  ${srcdir}/sorc/gdas.cd/modulefiles
#module load GDAS/hera


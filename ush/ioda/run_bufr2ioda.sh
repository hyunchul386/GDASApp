#!/bin/bash
#--------------------------------------------------------------------------------------
# run_bufr2ioda.sh
# This driver script will:
# - determine list of input BUFR files available
# - generate YAMLs from templates for each BUFR file
# - run BUFR2IODA.x and produce output IODA files
# usage:
#       run_bufr2ioda.sh YYYYMMDDHH /path/to/files.bufr_d/ /path/to/templates.yaml/ /path/to/output.ioda/
#
#--------------------------------------------------------------------------------------
if [[ $# -ne 5 ]] ; then
    echo "usage:"
    echo "      $0 YYYYMMDDHH gdas|gfs /path/to/files.bufr_d/ /path/to/templates.yaml/ /path/to/output.ioda/"
    exit 1
fi

# input parameters
CDATE=${CDATE:-$1}
RUN=${RUN:-$2}
BUFR_dir=${BUFR_dir:-$3}
YAML_template_dir=${YAML_template_dir:-$4}
out_dir=${out_dir:-$5}

# derived parameters
PDY=$(echo $CDATE | cut -c1-8)
cyc=$(echo $CDATE | cut -c9-10)

# get gdasapp root directory
readonly DIR_ROOT=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )/../.." && pwd -P)
BUFR2IODA=$DIR_ROOT/build/bin/bufr2ioda.x
BUFRYAMLGEN=$DIR_ROOT/ush/ioda/gen_bufr2ioda_yaml.py

# get list of BUFR files in input directory using glob
BUFR_files=$(ls $BUFR_dir/${RUN}.t${cyc}z.*.bufr_d)
if [ $? -ne 0 ]; then
    echo "No BUFR files found! in $BUFR_dir"
    exit 1
fi

# create output directory if it doesn't exist
mkdir -p $out_dir
if [ $? -ne 0 ]; then
    echo "cannot make $out_dir"
    exit 1
fi

# loop through available BUFR files
for f in $BUFR_files; do
    # get BUFR type from input BUFR file name
    BUFRbase=$(basename $f)
    BUFRtype_base=$(echo "${BUFRbase%.*.*}")
    BUFRtype=$(echo "${BUFRtype_base#*.*.}")
    echo "Now processing $BUFRtype"
    # get path to YAML template file for this BUFR file
    YAML_template=$YAML_template_dir/bufr_${BUFRtype}.yaml
    # check if the YAML template exists
    if [ ! -f $YAML_template ]; then
        echo "${YAML_template} does not exist, skipping ${BUFRtype}!"
        continue
    fi

    # input YAML file to the template parser
    cat > $out_dir/config_bufr_${BUFRtype}.yaml << EOF
obtype: $BUFRtype
input file: $f
output yaml file: $out_dir/bufr_${BUFRtype}.yaml
output dir: $out_dir
template yaml: $YAML_template
run: $RUN
PDY: $PDY
cyc: $cyc
EOF

    # now create YAML from the template
    $BUFRYAMLGEN --config $out_dir/config_bufr_${BUFRtype}.yaml

    # run BUFR2IODA for the created YAML file
    $BUFR2IODA $out_dir/bufr_${BUFRtype}.yaml

    # check if converter was successful
    if [ $? == 0 ]; then
      # remove YAMLs if success
      rm -rf $out_dir/config_bufr_${BUFRtype}.yaml
      rm -rf $out_dir/bufr_${BUFRtype}.yaml
    else
      # warn and keep YAMLs on failure
      echo "Problem running bufr2ioda.x for ${BUFRtype}"
      echo "See YAMLs in $out_dir"
      echo "  $out_dir/config_bufr_${BUFRtype}.yaml"
      echo "  $out_dir/bufr_${BUFRtype}.yaml"
    fi

done

# do extra stuff for the prepbufr file
#-----------
# first, check if the prepBUFR file exists
#if [ -f $BUFR_dir/${RUN}.t${cyc}z.prepbufr ]; then
#    echo "Will process $BUFR_dir/${RUN}.t${cyc}z.prepbufr"
#    # there are multiple YAMLs to run through
#    prepBUFR_YAMLS="prepbufr_adpupa prepbufr_adpsfc"
#fi

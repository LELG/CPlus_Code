#!/usr/bin/env bash

if [[ ! -z $PBS_SERVER ]] && [[ $PBS_SERVER == "bioinf-head.petermac.org.au" ]]; then
  export TERM=xterm
fi

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="config_file test_group param_set param_set_dir runs_per_ps"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=5
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of args provided."
fi

if [ ! -f $1 ]; then
  usage-exit "Could not find config file: "$1
else
  param_set_config_file=$1
fi

test_group="$2"
param_set="$3"
param_set_dir="$4"

# check that runs_per_ps is numeric
if $(echo $5 | grep -E -q '^[0-9]+$'); then
  runs_per_ps=$5
else
  usage-exit "runs_per_ps - numeric arg required, got "$2" instead"
fi

# check that simulation executable exists. No point doing any work if not.
if [ ! -f build/tumourSim ]; then
  echo "    ERR: tumourSim executable does not exist. Run build.sh and try again."
  exit 1
fi

# get appropriate padding for run_dir
# (this just finds the length of runs_per_ps as a string)
runpadding=${#runs_per_ps}

run_number=1
while [ $run_number -le $runs_per_ps ]; do
  # note need to 'pad' run directories for proper ordering (001, 002 ... etc)
  run_dir=$(printf "%s/%0*d" $param_set_dir $runpadding $run_number)

  if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
    # mkdir -p "$run_dir/data"
    # mkdir -p "$run_dir/plots"
  fi

  echo
  echo "------------------------------------------"
  echo
  echo "Starting simulation ..."
  echo
  echo "SIMULATION TEST GROUP: "$test_group
  echo "PARAMETER SET: "$param_set
  echo "RUN: "$run_number" of "$runs_per_ps
  echo
  echo "------------------------------------------"
  echo

  # run simulation with the full set of parameters
  build/tumourSim --run_number $run_number --run_dir $run_dir --config $param_set_config_file

  run_number=$((run_number+1))
done

#!/bin/bash

# if we're on the cluster, make sure correct PYTHONPATH is exported
if [[ ! -z $PBS_SERVER ]] && [[ $PBS_SERVER == "bioinf-head.petermac.org.au" ]]; then
  export PYTHONPATH=/usr/local/cluster/all_arch/python_libraries/production/lib/python2.7/site-packages
  export TERM=xterm
fi

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="config_file test_group param_set param_set_dir runs_per_param_set"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=5
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of args provided."
fi

# read parameters from config file
echo "================================================================================"
#printf "Reading parameters from config file ...\r"
param_set_config_file="$1"
test_group="$2"
param_set="$3"
param_set_dir="$4"
runs_per_param_set="$5"

# get appropriate padding for run_dir
# (this just finds the length of runs_per_param_set as a string)
runpadding=${#runs_per_param_set}

run_number=1
while [ $run_number -le $runs_per_param_set ]; do
  # note need to 'pad' run directories for proper ordering (001, 002 ... etc)
  run_dir=$(printf "%s/%0*d" $param_set_dir $runpadding $run_number)

  if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
    mkdir -p "$run_dir/data"
    mkdir -p "$run_dir/plots"
  fi

  echo "------------------------------------------"
  echo
  echo "Starting simulation ..."
  echo
  echo "SIMULATION TEST GROUP: "$test_group
  echo "PARAMETER SET: "$param_set #" of "$num_param_sets
  echo "RUN: "$run_number" of "$runs_per_param_set
  echo
  echo "------------------------------------------"
  echo

  echo "run_number = $run_number
run_dir = $run_dir" >> $param_set_config_file

  # run simulation with the full set of parameters
  build/testParser --config $param_set_config_file

  run_number=$((run_number+1))
done

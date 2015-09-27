#!/usr/bin/env bash

# Run treatment simulations on a local computer

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="num_runs"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=1
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of parameters provided"
fi

if $(echo $1 | grep -E -q '^[0-9]+$'); then
  num_runs=$1
else
  usage-exit "num_runs - expected numeric arg, got '"$1"' instead"
fi

# check that simulation executable exists. No point doing any work if not.
if [ ! -f build/tumourSimTreatment ]; then
  echo "    ERR: tumourSimTreatment executable does not exist. Run build.sh and try again."
  exit 1
fi

run_number=1
while [ $run_number -le $num_runs ]; do

  echo
  echo "------------------------------------------"
  echo
  echo "Starting TREATMENT simulation ..."
  echo
  echo "RUN: "$run_number" of "$num_runs
  echo
  echo "------------------------------------------"
  echo

  # run simulation
  mpirun -np 2 build/tumourSimTreatment > debugging_"$run_number".txt

  run_number=$((run_number+1))
done

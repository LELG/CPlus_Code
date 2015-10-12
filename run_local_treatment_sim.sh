#!/usr/bin/env bash

# Run treatment simulations on a local computer

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="num_runs num_proc"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=2
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of parameters provided"
fi

if $(echo $1 | grep -E -q '^[0-9]+$'); then
  num_runs=$1
else
  usage-exit "num_runs - expected numeric arg, got '"$1"' instead"
fi

if $(echo $2 | grep -E -q '^[0-9]+$'); then
  num_proc=$2
else
  usage-exit "num_runs - expected numeric arg, got '"$2"' instead"
fi

# check that simulation executable exists. No point doing any work if not.
if [ ! -f build/tumourSimTreatment ]; then
  echo "    ERR: tumourSimTreatment executable does not exist. Run build.sh and try again."
  exit 1
fi

test_group=$(date +"%a_%F_%H_%M_%S")

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
  mpirun -np $num_proc build/tumourSimTreatment --test_group $test_group --run_number $run_number

  run_number=$((run_number+1))
done

#!/usr/bin/env bash

# Run treatment simulations on a local computer

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="num_runs num_proc config_file"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=3
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

if [ ! -f $3 ]; then
  usage-exit "could not read config file: "$3
fi

# check that simulation executable exists. No point doing any work if not.
if [ ! -f build/tumourSimTreatment ]; then
  echo "    ERR: tumourSimTreatment executable does not exist. Run build.sh and try again."
  exit 1
fi

test_group=$(date +"%a_%F_%H_%M_%S")

# strip comments from config file
conf_file_no_comments=$(sed '/^#/'d $3)

# get param names from header row of config file,
# and save them to an array
IFS=','
param_names=($(echo "$file_no_comments" | awk 'NR==1'))
unset IFS

num_params=${#param_names[@]}

line=$(echo "$file_no_comments" | sed 1d )

config_file="$test_group.conf"

# save this set of parameter values to an array
IFS=','
vals=($(echo "$line"))
unset IFS

printf "==> Writing parameter set to config file ...\r"

i=0
while [ $i -lt $num_params ]; do
  param_name=${param_names[i]}
  param_val="${vals[i]}"

  if [ "$param_val" = "false" ] || [ -z "$param_val" ]; then
    # param val is "false"/empty string
    let i++
    continue
  fi

  # check if param is just a boolean switch
  if [ "$param_val" = "true" ]; then
    echo "$param_name = true" >> $config_file
  else
    # check for multi-valued param value
    if [ $(wc -w <<< "$param_val") -gt 1 ]; then
      # split multiple values into array
      subvals=($(echo "$param_val"))
      # argparse requires the values on separate lines
      for subval in "${subvals[@]}"; do
        echo "$param_name = $subval" >> $config_file
      done
    else
      # single-valued param val
      echo "$param_name = $param_val" >> $config_file
    fi
  fi

  let i++
done

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
  mpirun -np $num_proc build/tumourSimTreatment --test_group $test_group --run_number $run_number --config $config_file

  run_number=$((run_number+1))
done

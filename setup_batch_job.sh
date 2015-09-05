#!/usr/bin/env bash

# Generate script for running batch (e.g. PBS) jobs
# and create some directories for storing results

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="[-m MAIL_ADDRESS] [-w WALLTIME] [-q QUEUE_NAME] testgroup_name runs_per_ps config_file"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

# process flag options using getopts

OPTIND=1    # Reset in case getopts has been used previously in the shell.

# initialise batch variables
walltime="1:00:00"
queue=""
mailopt=""
mailadd=""

while getopts "h?w:m:q:" opt; do
    case "$opt" in
    h|\?)
        usage-exit
        ;;
    w)  walltime=$OPTARG
        ;;
    m)  mailopt="#PBS -m ae"
	mailadd="#PBS -M "$OPTARG
	;;
    q)  queue="#PBS -q "$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# now check for positional arguments

num_expected_args=3
if [ $# -ne $num_expected_args ]; then
  usage-exit
fi

# check that supplied testname will make a valid directory name
# i.e., that it does not contain slash or backslash
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
  usage-exit "Test group name cannot contain slash/backslash"
else
  test_group=$1
fi

# check that runs_per_ps is numeric
if $(echo $2 | grep -E -q '^[0-9]+$'); then
  runs_per_ps=$2
else
  usage-exit "runs_per_ps - numeric arg required, got "$2" instead"
fi

if [ ! -f $3 ]; then
  usage-exit "Could not find config file: "$3
else
  tg_config_file=$3
fi


# make directory for today's date, unless it already exists
today=$(date +'%Y-%m-%d')
today_dir="results/$today"
if [ ! -d "$today_dir" ]
then
  printf "=== Creating results directory for "$today" ...\r"
  mkdir -p "$today_dir"
  printf "=== Creating results directory for "$today" ... done.\n"
fi

# make main directory for this test group
test_group_dir="$today_dir/$test_group"
if [ ! -d $test_group_dir ]
then
  printf "=== Creating test group directory: "$test_group_dir" ...\r"
  mkdir -p $test_group_dir
  printf "=== Creating test group directory: "$test_group_dir" ... done.\n"
else
  echo "    WARNING: results for test group "$test_group" already exist"
  i=1
  while [ -d "$test_group_dir($i)" ] ; do
    let i++
  done
  test_group_dir="$test_group_dir($i)"
  test_group="$test_group($i)"
  printf "=== Creating test group directory: "$test_group_dir" ...\r"
  mkdir -p $test_group_dir
  printf "=== Creating test group directory: "$test_group_dir" ... done.\n"
fi

if [ ! -f $test_group_dir/$tg_config_file ]; then
  printf "=== Copying config file to $test_group_dir ...\r"
  cp $tg_config_file $test_group_dir
  printf "=== Copying config file to $test_group_dir ... done.\n"
fi

# strip comments from config file
file_no_comments=$(sed '/^#/'d $tg_config_file)
# strip header
param_sets=$(echo "$file_no_comments" | sed 1d)

# get number of non-header lines in param file
# xargs trims whitespace from return value of wc
num_param_sets=$(echo "$param_sets" | wc -l | xargs)

# get param names from header row of config file,
# and save them to an array
IFS=','
param_names=($(echo "$file_no_comments" | awk 'NR==1'))
unset IFS

num_params=${#param_names[@]}
param_set_padding=${#num_param_sets}

# start constructing PBS script
pbs_script="PBS-$test_group.sh"
echo "=== Creating PBS script for \`$test_group\` ..."

# first, send necessary variables
cat >> $pbs_script << _endmsg
#!/bin/bash
#PBS -N $test_group
#PBS -l walltime=$walltime
#PBS -o $test_group.log
#PBS -j oe
#PBS -t 1-$num_param_sets
$queue
$mailopt
$mailadd

this_script="$pbs_script"
tg_config_file="$tg_config_file"
param_sets="$param_sets"
param_names=(${param_names[@]})
num_params=${#param_names[@]}
num_param_sets=$num_param_sets
param_set_padding=${#num_param_sets}
runs_per_ps=$runs_per_ps
test_group="$test_group"
test_group_dir="$test_group_dir"
_endmsg

# then send the body of the script, verbatim
cat >> $pbs_script << '_endmsg'

cd $PBS_O_WORKDIR

# use PBS array id to iterate over param sets
param_set=$PBS_ARRAYID

# get correct param set
line=$(echo "$param_sets" | awk -v linenum=$param_set 'NR==linenum')

# save this set of parameter values to an array
IFS=','
vals=($(echo "$line"))
unset IFS

# pad directories to ensure proper ordering
param_set_dir=$(printf "%s/%0*d" $test_group_dir $param_set_padding $param_set)

if [ ! -d $param_set_dir ]; then
  printf "=== Creating param set directory: "$param_set_dir" ...\r"
  mkdir -p $param_set_dir
  printf "=== Creating param set directory: "$param_set_dir" ... done.\n"
fi

# each parameter set will get stored to its own file
# (in a boost::program_options compatible format)
param_set_config_file="$param_set_dir/$test_group-$param_set.cpp.conf"

printf "=== Writing parameter set to config file ...\r"

i=0
while [ $i -lt $num_params ]; do
  param_name=${param_names[i]}
  param_val="${vals[i]}"

  if [ "$param_name" = "scale" ] || [ "$param_name" = "mscale" ]; then
    # these parameters are not required in C++ simulation
    let i++
    continue
  fi

  if [ "$param_val" = "false" ] || [ -z "$param_val" ]; then
    # param val is "false"/empty string
    let i++
    continue
  fi

  # check if param is just a boolean switch
  if [ "$param_val" = "true" ]; then
    echo "$param_name = true" >> $param_set_config_file
  else
    # check for multi-valued param value
    if [ $(wc -w <<< "$param_val") -gt 1 ]; then
      # split multiple values into array
      subvals=($(echo "$param_val"))
      # argparse requires the values on separate lines
      for subval in "${subvals[@]}"; do
        echo "$param_name = $subval" >> $param_set_config_file
      done
    else
      # single-valued param val
      echo "$param_name = $param_val" >> $param_set_config_file
    fi
  fi

  let i++
done

# echo misc variables to the file
echo "test_group = $test_group
param_set = $param_set
test_group_dir = '$test_group_dir'
param_set_dir = '$param_set_dir'" >> $param_set_config_file

printf "=== Writing parameter set to config file ... done.\n"

./run_param_set.sh $param_set_config_file $test_group $param_set $param_set_dir $runs_per_ps

# when run, this script copies itself to the test group directory
cp "$this_script" "$test_group_dir"
_endmsg

echo "=== PBS script created as $pbs_script."

#!/usr/bin/env bash

# Generate script for running batch (e.g. PBS) jobs
# for the treatment version of simulation

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="[-m MAIL_ADDRESS] [-w WALLTIME] [-q QUEUE_NAME] num_runs num_proc"

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

num_expected_args=2
if [ $# -ne $num_expected_args ]; then
  usage-exit
fi

# check that num_runs is numeric
if $(echo $1 | grep -E -q '^[0-9]+$'); then
  num_runs=$1
else
  usage-exit "num_runs - expected numeric arg, got '"$1"' instead"
fi

# check that num_proc is numeric
if $(echo $2 | grep -E -q '^[0-9]+$'); then
  num_proc=$2
else
  usage-exit "num_proc - expected numeric arg, got '"$2"' instead"
fi

# start constructing PBS script
test_group=$(date +"%a_%F_%H_%M_%S")
pbs_script="PBS-treatment-"$test_group".sh"
echo "==> Creating PBS script $pbs_script ..."

# first, send necessary variables
cat >> $pbs_script << _endmsg
#!/bin/bash
#PBS -N treatment_$test_group
#PBS -l walltime=$walltime
#PBS -l nodes=1:ppn=$num_proc
#PBS -o treatment_$test_group.log
#PBS -j oe
#PBS -t 1-$num_runs
$queue
$mailopt
$mailadd

this_script="$pbs_script"
num_runs=$num_runs
test_group="$test_group"
_endmsg

# then send the body of the script, verbatim
cat >> $pbs_script << '_endmsg'

cd $PBS_O_WORKDIR

# use PBS array id to iterate over runs
run_number=$PBS_ARRAYID


# check that simulation executable exists. No point doing any work if not.
if [ ! -f build/tumourSimTreatment ]; then
  echo "    ERR: tumourSimTreatment executable does not exist. Run build.sh and try again."
  exit 1
fi

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
mpirun build/tumourSimTreatment --test_group $test_group --run_number $run_number
_endmsg

echo "==> PBS script created"

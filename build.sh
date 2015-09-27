#!/usr/bin/env bash

script_parameters="[--treatment]"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

if [ $# -ne 0 -a $# -ne 1 ]
then
  usage-exit "Incorrect number of parameters provided"
elif [ $# -eq 0 ]
then
  treatment=false
elif [ "$1" = --treatment ]
then
  treatment=true
else
  usage-exit "Unknown option '"$1"'"
fi

if [ "$treatment" = true ]
then
  echo "==> Starting simulation build (version with TREATMENT)"
else
  echo "==> Starting simulation build"
fi

mkdir build 2> /dev/null

if [ $? -eq 0 ]
then
  echo "==> Created build directory"
else
  echo "ERR: Build directory already exists. Run \`clean.sh\`. Aborting." >&2
  exit 1
fi

echo "==> Running CMake"

touch build/build.log

if [ $? -ne 0 ]; then
  echo "ERR: could not create build log file. Aborting." >&2
  exit 1
fi

( cd build && CC=gcc CXX=g++ cmake .. &> build.log )

if [ $? -ne 0 ]; then
  echo "ERR: CMake failed. See build/build.log for details." >&2
  exit 1
fi

echo "==> Running CMake --build ..."

printf "\nbuilding CMake target\n\n" >> build/build.log
if [ "$treatment" = true ]
then
  ( cd build && make tumourSimTreatment >> build.log 2>&1 )
else
  ( cd build && make tumourSim >> build.log 2>&1 )
fi

if [ $? -ne 0 ]; then
  echo "ERR: CMake build failed. See build/build.log for details." >&2
  exit 1
fi

echo "==> Simulation built. See README.md for intructions on running."
exit 0

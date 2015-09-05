#!/usr/bin/env bash

echo "=== Starting simulation build"

mkdir build 2> /dev/null

if [ $? -eq 0 ]
then
  echo "=== Created build directory"
else
  echo "ERR: Build directory already exists. Run \`clean.sh\`. Aborting." >&2
  exit 1
fi

echo "=== Running CMake"

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

echo "=== Running CMake --build ..."

cmake --build build >>build/build.log 2>&1

if [ $? -ne 0 ]; then
  echo "ERR: CMake build failed. See build/build.log for details." >&2
  exit 1
fi

echo "=== Simulation built. See README.md for intructions on running."
exit 0

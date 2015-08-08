#!/bin/bash

if [ -d build ]; then
  # remove build directory
  rm -rf build
  echo "=== Build directory removed"
  exit 0
else
  echo "=== Nothing to clean up"
  exit 1
fi

#!/bin/bash

womtool_jar="configs/womtool-88.jar"

# check_miniwdl() {
#   local wdl=$1
#   echo "> Running miniwdl check on $wdl:"
#   miniwdl check --no-quant-check --strict "$wdl"
#   echo ""
# }

check_womtool() {
  local wdl=$1
  echo "> Running womtool validate on $wdl:"
  java -jar "$womtool_jar" validate "$wdl"
  echo ""
}

# Loop over all main WDL files.
for wdl in *.wdl; do
  # check_miniwdl "$wdl"
  check_womtool "$wdl"
done
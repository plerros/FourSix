#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"method": "ant",//g' {} \;
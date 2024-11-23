#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

find "$SCRIPT_DIR/s_anneal" -type  f -name "*.json" -exec sed -i 's/{"L": 200, "alpha":5.0, "beta":0.2}/{"L": 400, "alpha":5.0, "beta":0.2}/g' {} \;
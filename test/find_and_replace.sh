#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

find "$SCRIPT_DIR/ant" -type  f -name "*.json" -exec sed -i 's/{"L": 10, "alpha":5.0, "beta":0.2, "xi": 1.0, "psi": 3.0, "lambda": 0.5, "kappa": 10}/{"L": 100, "alpha":5.0, "beta":0.2, "xi": 1.0, "psi": 3.0, "lambda": 0.5, "kappa": 10}/g' {} \;
#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/{"L":.*[0-9]*, "alpha":.*[0-9]*.[0-9]*, "beta":.*[0-9]*.[0-9]*}/{"L": 100, "alpha": 5.0, "beta": 0.2}/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"L":[0-9]*}/"L": 200}/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"L": [0-9]*}/"L": 200}/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"L":[0-9]*,/"L": 200,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"L": [0-9]*,/"L": 200,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"alpha":[0-9]*.[0-9]*,/"alpha": 5.0,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"alpha": [0-9]*.[0-9]*,/"alpha": 5.0,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"beta":[0-9]*.[0-9]*}/"beta": 0.2}/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"beta": [0-9]*.[0-9]*}/"beta": 0.2}/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"beta":[0-9]*.[0-9]*,/"beta": 0.2,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"beta": [0-9]*.[0-9]*,/"beta": 0.2,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"xi":[0-9]*.[0-9]*,/"xi": 1.0,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"xi": [0-9]*.[0-9]*,/"xi": 1.0,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"psi":[0-9]*.[0-9]*,/"psi": 3.0,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"psi": [0-9]*.[0-9]*,/"psi": 3.0,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"lambda":[0-9]*.[0-9]*,/"lambda": 0.5,/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"lambda": [0-9]*.[0-9]*,/"lambda": 0.5,/g' {} \;

find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"kappa":[0-9]*}/"kappa": 10}/g' {} \;
find "$SCRIPT_DIR/autodetect" -type  f -name "*.json" -exec sed -i 's/"kappa": [0-9]*}/"kappa": 10}/g' {} \;
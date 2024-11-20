#!/bin/bash

# detect if it falls in a very long loop
hyperfine -m 200 './FourSix -i ../test/example_instances_rev1/example_instances/cgshop2025_examples_ortho_20_b099d1fe.instance.json' --show-output

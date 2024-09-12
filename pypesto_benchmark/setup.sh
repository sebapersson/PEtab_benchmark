#!/usr/bin/env bash
python3 -m virtualenv venv --clear
source ./venv/bin/activate
pip install -r requirements.txt
git clone --depth 1 https://github.com/benchmarking-initiative/Benchmark-Models-PEtab.git

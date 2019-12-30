#!/usr/bin/env bash

cleanCompss;

# commented commands
# -g -t \

runcompss \
    --lang=python \
    --python_interpreter=python3 \
    --pythonpath=. \
    ./main.py
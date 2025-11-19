#!/bin/bash 

params=$1 
results=$2 

export PYTHONPATH=/Users/swyant/local_software/Dakota/share/dakota/Python 
uv run python DakotaSubDriver.py $params $results


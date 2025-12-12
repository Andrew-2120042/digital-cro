#!/bin/bash
cd "/Users/nareshnallabothula/digital cro/drug-discovery"
export PATH=~/bin:$PATH
python test_quick.py > smoke_output.txt 2>&1
echo "Exit code: $?" >> smoke_output.txt

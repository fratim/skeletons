#!/bin/bash
set -e

echo "Execution started"

# Execute step one
for bz in {0..1}
do
  for by in {0..1}
  do
    for bx in {0..1}
      do
        python scripts/connectome.py $bz $by $bx
    done
  done
done

echo "Step 1 finished"

#!/bin/bash

# Get the list of scripts, sorted by the numerical order
scripts=$(ls *.py *.R 2>/dev/null | grep -E '^[0-9]' | sort -V)

# Check if the scripts variable is empty
if [ -z "$scripts" ]; then
    echo "No Python or R scripts found."
    exit 1
fi

# Loop through each script and execute it based on its extension
for script in $scripts; do
    case "$script" in
        *.py)
            echo "Running Python script: $script"
            python3 "$script"
            ;;
        *.R)
            echo "Running R script: $script"
            Rscript "$script"
            ;;
        *)
            echo "Unknown file type: $script"
            ;;
    esac
done

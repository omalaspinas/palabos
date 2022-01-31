#!/bin/bash

if [ "$1" == "-n" ]; then
    echo "Running dry run. Warnings are errors."
    find examples src utility coupledSimulators -type d \( -path '*build*' \) -prune -false \
        -o -name '*.cpp' -o -name '*.h' -o -name '*.hh' | xargs clang-format -n -Werror

elif [ "$1" == "-i" ]; then
# We only include src, examples directories
# We don't want to add external libs
# Also we remove the build directories
# Finally we only want to reformat .h, .hh and .cpp files
# Clang format is run with -i option to change files in place
    echo "Formatting in place."
    find examples src utility coupledSimulators -type d \( -path '*build*' \) -prune -false \
        -o -name '*.cpp' -o -name '*.h' -o -name '*.hh' | xargs clang-format -i
else
    echo "An option is mandatory: can be -n for dry run or -i for in place."
fi

#!/bin/bash

# Check if the argument (test number) is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <test_number>"
  exit 1
fi

test_number="$1"

# Change to the 'src' directory
cd src || exit 1

# Compile the code
make || exit 1

# Copy the executable to the 'checker' directory
cp tema1_par ../checker || exit 1

# Change to the 'checker' directory
cd ../checker || exit 1

# Run the test with the provided test number
./tema1_par inputs/in_"$test_number".ppm ../outputs/out"$test_number".ppm 2 || exit 1

# Success
echo "Test $test_number completed successfully"
exit 0

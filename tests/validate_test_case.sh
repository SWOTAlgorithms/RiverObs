#!/usr/bin/env bash

# Run the python test suite
pytest -v test_riverobs.py

rm -rf ncdiffs
mkdir -p ncdiffs

if test -f "test_case/rt.nc"; then
  ./ncdiff_netcdfs.sh reference_case/rt.nc test_case/rt.nc ncdiffs/rt/
  ./ncdiff_netcdfs.sh reference_case/pv.nc test_case/pv.nc ncdiffs/pv/
else
  echo "Test case missing. Must run test case before validating it."
fi
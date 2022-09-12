#!/usr/bin/env bash

# Run the python test suite
pytest -v test_riverobs.py

rm -rf ncdiffs
mkdir -p ncdiffs
./ncdiff_netcdfs.sh reference_case/rt.nc test_case/rt.nc ncdiffs/rt/
./ncdiff_netcdfs.sh reference_case/pv.nc test_case/pv.nc ncdiffs/pv/

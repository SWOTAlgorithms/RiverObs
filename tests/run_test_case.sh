#!/usr/bin/env bash

PIXC_FILE=inputs/SWOT_L2_HR_PIXC_001_578_086L_20190221T145823_20190221T145827_DG00_01.nc

AUX_FILE=inputs/SWOT_Param_L2_HR_RiverTile_20200101T000000_21000101T000000_20220810T000000_v201.rdf


mkdir -p test_case
$SWOT_SDS_BASE/RiverObs/src/bin/swot_pixc2rivertile.py \
    $PIXC_FILE test_case/rt.nc test_case/pv.nc $AUX_FILE


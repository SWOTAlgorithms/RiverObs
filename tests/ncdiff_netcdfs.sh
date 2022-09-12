#!/bin/bash
# Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
# Government sponsorship acknowledged.
# All rights reserved.
#
# Author(s): Bryan Stiles

if [ "$#" -lt 3 ]; then
echo "Usage ncdiff_netcdfs.sh ncfile1 ncfile2 output_dir [normalize]"
echo " If normalize==1 then normalize \both files so that each variable in the reference file is zero mean and RMS=1"
echo " Infile 2 is presumed to be the reference file"
exit
fi


infile1=$1
infile2=$2
outdir=$3

norm=0
if [ "$#" -eq 4 ]; then
norm=$4
fi

echo "Comparing $infile1 to $infile2 and placing results in $outdir directory"

# make temporary files without nans and with _FillValue=1E36
echo "making temporary files without nans and with _FillValue=1E36"
tmp1=infile1_nonan.nc
tmp2=infile2_nonan.nc
mkdir $outdir
cp $infile1 $outdir/$tmp1
cp $infile2 $outdir/$tmp2
cd $outdir

ncatted -a _FillValue,,o,f,NaN $tmp1
ncatted -a _FillValue,,m,f,1.0e36 $tmp1

ncatted -a _FillValue,,o,f,NaN $tmp2
ncatted -a _FillValue,,m,f,1.0e36 $tmp2

# normalize files if desired
if [ "$norm" -eq 1 ]; then
ncwa $tmp1 mean1.nc
ncwa $tmp2 mean2.nc
ncbo --op_typ='-'  $tmp1 mean2.nc anomaly1.nc
ncbo --op_typ='-'  $tmp2 mean2.nc anomaly2.nc
ncwa -y rmssdn anomaly1.nc std1.nc
ncwa -y rmssdn anomaly2.nc std2.nc
ncbo --op_typ='/' anomaly1.nc std2.nc infile1_norm.nc
ncbo --op_typ='/' anomaly2.nc std2.nc infile2_norm.nc
tmp1=infile1_norm.nc
tmp2=infile2_norm.nc
fi
# make difference files
echo "making difference files"

ncdiff $tmp1 $tmp2 diff.nc
ncdiff $tmp2 $tmp2 perfectdiff.nc


# make mean diff files and text report
echo "making mean difference files and report"
ncwa diff.nc mean_diff.nc
ncwa perfectdiff.nc mean_perfectdiff.nc

ncdump mean_diff.nc > tmpdump
ncdump mean_perfectdiff.nc > tmp2dump
substring.csh group GROUP tmpdump
diff tmpdump tmp2dump > tmpdump3
substring.csh '<' "" tmpdump3
grep -v '>' tmpdump3 > mean_diff.txt

# make mean square diff files and text report
echo "making rms difference files and report"
ncwa -y rms diff.nc rms_diff.nc
ncwa -y rms perfectdiff.nc rms_perfectdiff.nc
ncdump rms_diff.nc > tmpdump
ncdump rms_perfectdiff.nc > tmp2dump
substring.csh group GROUP tmpdump
diff tmpdump tmp2dump > tmpdump3
substring.csh '<' "" tmpdump3
grep -v '>' tmpdump3 >  rms_diff.txt
\rm tmpdump tmp2dump tmpdump3

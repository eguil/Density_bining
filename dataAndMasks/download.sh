#!/bin/tcsh

# File written to download all OAFlux_V3*.nc monthly files
# PJD 22 Jun 2014	- Adapted from ERSSTV3b for OAFlux_v3 files

set date=`date +%y%m%d`
set wodpath=/work/durack1/Shared/obs_data/WOD13
cd ${wodpath}
\rm -r -f ${date} ; # Purge if exists
\mkdir ${date}
cd ${date}
# Download files
\wget -o ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/basinmask_01.msk
\wget -a ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/basinmask_04.msk
\wget -a ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/landsea_01.msk
\wget -a ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/landsea_04.msk
\wget -a ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/mixnumber_01.msk
\wget -a ../${date}_log.txt -nv -nH --cut-dirs=3 http://data.nodc.noaa.gov/woa/WOA13/MASKS/mixnumber_04.msk
# Gunzip all files into a temp space: ls -1 *nc.gz | sed "s/\(.*\)\.gz/gunzip -c \1.gz > \1/" | sh
# Process individual files
cd ${wodpath}
echo "Starting python netcdf conversion.."
#uvcdat6 make_masks.py

: << --
# Set symlinks
#echo "Create symlink to new data.."
#\ls -1 ${date}_*.nc | sed "s/^[^_]*_\(.*\)_.*\.nc/ln -fs '&' '\1.nc'/" | sh ; # local ln creation, truncating YYYYMM-YYYYMM from filename
# Zip up source data
echo "Archive source data and purge temp ${date} directory.."
\rm -f ${date}.tar.bz2
\tar -cjf ${date}.tar.bz2 ${date} ; # Archive using bzip2 compression
# Extract using >tar -xjf ${date}.tar.bz2
# Grab descriptor files and purge directory
mv -f ${date}/OAFlux_TechReport_3rd_release.pdf .
mv -f ${date}/readme-monthly.txt .

#rm -rf ${date}
echo "${date}_OAFlux_v3.nc download complete.."
# Conditionally clean up existing versions of files and archive
if ( $1 != "" ) then
	\tar -cjf ${1}_archive.tar.bz2 ${1}*.* ; # Archive using bzip2 compression
	\rm -r -f ${1}.tar.bz2
	\rm -r -f ${1}_OAFlux_V3_*.nc
	\rm -r -f ${1}_log.txt
endif
--
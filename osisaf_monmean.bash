#!/usr/bin/env bash 
echo "Works, but not needed. Monthly means can be downloaded directly from OSI SAF FTP: osisaf.met.no:/reprocessed/ice/conc/v3p0/monthly/"
exit 1
module load cdo
module load nco
INDIR="/albedo/work/user/vludwig/05_SINXS/02_DATA/00_MISC/osisaf/"
OUTDIR="/albedo/work/user/vludwig/05_SINXS/02_DATA/00_MISC/osisaf/monthly/"

for YEAR in {1999..2023}; do
	for MONTH in 01 02 03 04 05 06 07 08 09 10 11 12; do
			if [ -d $INDIR/$YEAR/$MONTH ]; then 
				echo "########################"
				echo "### STARTING $YEAR/$MONTH ###"
				echo "########################"
			else
				echo "No files for $INDIR/$YEAR/$MONTH; continue "
				continue
			fi
			if [ "$YEAR" -lt 2021 ]; then
				OUTFILE=ice_conc_nh_ease2-250_cdr-v3p0_$YEAR$MONTH.nc
				cdo timmean -cat $INDIR/$YEAR/$MONTH/'ice_conc_nh_ease2-250_cdr-v3p0_*.nc'  $OUTDIR/MONMEAN.nc
				cdo selname,ice_conc $OUTDIR/MONMEAN.nc $OUTDIR/$OUTFILE
			else
				OUTFILE=ice_conc_nh_ease2-250_icdr-v3p0_$YEAR$MONTH.nc
				cdo timmean -cat $INDIR/$YEAR/$MONTH/'ice_conc_nh_ease2-250_icdr-v3p0_*.nc'  $OUTDIR/MONMEAN.nc
				cdo selname,ice_conc $OUTDIR/MONMEAN.nc $OUTDIR/$OUTFILE
			fi
			rm $OUTDIR/MONMEAN.nc
			ncrename -v ice_conc,sea_ice_concentration $OUTDIR/$OUTFILE
	done
done

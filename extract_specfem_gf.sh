#!/usr/bin/env bash

specfem_dir=$1
station=$2
channel=$3
temp_output_dir=$4

# copy output files (only txt output)
mkdir -p $station
cp $specfem_dir/output_* $station/

# concatenate specfem output (typically one file per rank > one big file)
cat $specfem_dir/all_seismograms_*.bin > $station/$station.bin

# extract, interpolate to lower sampling rate
# depending on the shell this is used on, may have to include additional
# commands to activate python environment like source <env>/bin/activate for virtual env
source ~/code/miniconda3/etc/profile.d/conda.sh
conda activate noisi
mpirun -np 4 python parse_specfem_output.py $station/$station.bin

# concatenate again
cat $temp_output_dir/$station..$channel/*.bin_* > $station/$station..$channel.bin

# finally, transform from unformatted binary to hdf5
# best to do this separately in a task farm (to save resources)
# python .bin_to_h5.py $station/$station..$channel.bin ../sourcegrid_densergrid.npy

rm -rf $temp_output_dir



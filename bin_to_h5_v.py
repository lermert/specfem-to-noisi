# ToDo check whether these are all the metadata we sensibly need.
import h5py
import os
import sys
import numpy as np

# User input: -------------------------------------------
# -------------------------------------------------------
nbytes_stationname = 512
size_of_float = 4  # which precision was used by specfem
data_quantity = 'VEL'  # 'DIS','VEL','ACC'
source_component = "z"
# -------------------------------------------------------
# -------------------------------------------------------
f_in = sys.argv[1]
f_sources = sys.argv[2]
nbytes_total = os.path.getsize(f_in)

f_out_name = os.path.splitext(f_in)[0] + '.h5'
f_sources = np.load(f_sources)
f_out = h5py.File(f_out_name, "w")
input_file = open(f_in, 'rb')
# Get metadata
input_file.seek(nbytes_stationname)
ntimesteps = np.fromfile(input_file,
                         dtype='f' + str(size_of_float), count=1)[0]
ntimesteps = int(ntimesteps)
Fs = round(np.fromfile(input_file,
                       dtype='f' + str(size_of_float), count=1)[0], 6)
# Record lengths: station name plus two header values plus length of data array
nbytes_trace = nbytes_stationname + (2 + ntimesteps) * size_of_float
# Number of records actually contained
ntraces = int(nbytes_total / nbytes_trace)
if nbytes_total % nbytes_trace != 0:
    raise ValueError("Error decoding bin format. Check size of float and\
bytes per station name.")
print('This file contains %g Traces.' % ntraces)

# Reference station:
refstation = os.path.basename(sys.argv[1])
refstation = os.path.splitext(refstation)[0]

# DATASET NR 1: STATS
if 'stats' not in list(f_out.keys()):
    stats = f_out.create_dataset('stats', data=(0,))
    stats.attrs['reference_station'] = refstation
    stats.attrs['data_quantity'] = data_quantity
    stats.attrs['ntraces'] = ntraces
    stats.attrs['Fs'] = Fs
    stats.attrs['nt'] = int(ntimesteps)
    stats.attrs['fdomain'] = False

# DATASET NR 2: Source grid
if 'sourcegrid' not in list(f_out.keys()):
    sources = f_out.create_dataset('sourcegrid', data=f_sources[0:2])

# DATASET Nr 3: Seismograms itself
traces = f_out.create_dataset('data',
                              (ntraces, ntimesteps), dtype=np.float32)

# jump to the beginning of the trace in the binary file
input_file.seek(0)
i = 0
print('Starting to read seismograms from: %s' % sys.argv[1])
while i < ntraces:
    if i % 10000 == 0:
        print('Converted %g of %g traces' % (i, ntraces))
    # read station name, copy to output file
    staname = input_file.read(nbytes_stationname)
    staname = str(staname.decode('utf-8')).strip()
    print(staname)
    # These are only read to jump over the entries
    nt_temp = int(np.fromfile(input_file, dtype='f' + str(size_of_float),
                              count=1)[0])
    Fs_temp = np.fromfile(input_file, dtype='f' + str(size_of_float), count=1)[0]
    # Get the index of that station
    # This links it with the right source coordinate pair.
    staindex = int(staname.split('.')[1])
    values = np.fromfile(input_file, dtype='f' + str(size_of_float),
                         count=ntimesteps)
    # Save in traces array
    traces[staindex, :] += values
    i += 1
f_out.close()

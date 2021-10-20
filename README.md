# specfem-to-noisi
Some scripts to prepare input for the global version of specfem, and parse the output into noisi-friendly hdf5 files.
[noisi](https://github.com/lermert/noisi) is a tool to calculate noise cross-correlations and source sensitivity kernels


Dependencies: scipy, numpy, mpi4py, h5py, obspy. (noisi for generating source grid file and plotting).
To use simulations from specfem3d_globe for noise correlations:
- prepare two files (see example directory for examples of these two types of file)
	- a stations list containing comma separated values of name, lat, lon, elev_m, sensor_depth_m
	- a noise source grid (sourcegrid.npy) file using noisi

- use the scripts in prepare_sem_input.py to prepare the specfem input: 
	- One FORCESOLUTION file per station
	- STATIONS file where each STATIONS is a source grid point

- edit the specfem Par_file (an example is in the example directory). Note that these scripts assume we will use binary output files in specfem3d_globe. This is the simplest option, because they can be written very fast by specfem.
- run specfem3d_globe using the STATIONS list, Par_file and FORCESOLUTION as inputs
- edit parse_specfem_output.py using the output_solver.txt file writen by specfem.
- concatenate the specfem output, and parse it using extract_specfem_gf.sh
- convert the binary format to hdf5 using bin_to_h5_v.py

If all these steps worked well, then plotting the wavefield with python -c 'from noisi import WaveField; wf = WaveField("SSB/SSB..MXZ.h5"); wf.plot_snapshot(180)' should look similar to the example plot in the example folder. :)

Doesn't work? Feel free to submit an issue or otherwise contact me.
import numpy as np
from pandas import read_csv
import os


hdur_pointsource = 2.0
outdir = 'specfem_input'
# force source directions:
components = ['Z']
sourcegridfile = "example/sourcegrid.npy"
# station list:
stationcsvfile = "example/stationlist.csv"
# convert station latitude or not? IRIS latitudes: Do Not convert
convert_stala = False
# use topography?
use_topo = 1


def wgs84():

    # semi-major axis, in m
    a = 6378137.0

    # semi-minor axis, in m
    b = 6356752.314245

    # inverse flattening f
    f = a / (a - b)

    # squared eccentricity e
    e_2 = (a ** 2 - b ** 2) / a ** 2

    return(a, b, e_2, f)


def geocent_to_geograph(theta):
    # convert geocentric latitudes to geographic assuming WGS84
    # the other way around
    e2 = wgs84()[2]
    theta = np.rad2deg(np.arctan(np.tan(np.deg2rad(theta)) / (1 - e2)))
    return(theta)


def grid_to_stations_file(outdir, grid):
    """
    Write noisesource grid as STATIONS list that specfem and axisem3d can read
    """

    fid = open(os.path.join(outdir, 'STATIONS'), 'w')
    for i in range(len(grid[0, :])):
        fid.write('%08g SRC %10.8f  %10.8f 0.0 0.0\n'
                  % (i, grid[1, i], grid[0, i]))

    fid.close()


def stations_to_forcesolutions(stationlist, hdur, outdir,
                               channel='MXZ',
                               use_topo=False, use_buried=False,
                               source_type='Gauss', factor_force=1.e10):

    stationlist = read_csv(stationlist)
    os.makedirs(outdir, exist_ok=1)

    for i in range(len(stationlist)):
        station = stationlist.iloc[i].sta
        network = stationlist.iloc[i].net
        if convert_stala:
            latitude = geocent_to_geograph(stationlist.iloc[i].lat)
        else:
            latitude = stationlist.iloc[i].lat
        longitude = stationlist.iloc[i].lon
        station_id = network + '.' + station + '..' + channel
        print(station_id)
        if use_topo:
            elevation = stationlist.iloc[i].elev_m / 1000.
        else:
            elevation = 0.
        if use_buried:
            depth = stationlist.iloc[i].sensor_depth_m / 1000.
            depth = depth - elevation
        else:
            if not use_topo:
                depth = 0.
            else:
                depth = -elevation

        if source_type == 'Gauss':
            stftype = 0
        elif source_type == 'Ricker':
            stftype = 1

        for component in components:
            eventfid = open(os.path.join(outdir,
                                         'FORCESOLUTION' + '_' + station +
                                         '_' + component), 'w')
            eventfid.write('FORCE 001 \n')
            eventfid.write('time shift:    0.0000   \n')
            eventfid.write('half duration:    %s   \n' % str(hdur))
            eventfid.write('latitude:    %s   \n' % str(latitude))
            eventfid.write('longitude:    %s   \n' % str(longitude))
            eventfid.write('depth:    %s   \n' % str(depth))
            eventfid.write('source time function:  %g \n' % stftype)
            if component == 'Z':
                f = [0.0, 0.0, -1.0]
            elif component == 'E':
                f = [1.0, 0.0, 0.0]
            elif component == 'N':
                f = [0.0, 1.0, 0.0]
            eventfid.write('factor force source:       %g \n' % factor_force)
            eventfid.write('component dir vect source E:   %g   \n' % f[0])
            eventfid.write('component dir vect source N:   %g   \n' % f[1])
            eventfid.write('component dir vect source Z_UP:   %g   \n' % f[2])


if __name__ == '__main__':

    grid = np.load(sourcegridfile)
    stations = os.path.join(stationcsvfile)
    stations_to_forcesolutions(stations, hdur_pointsource, outdir,
                               use_buried=1, use_topo=use_topo)
    grid_to_stations_file(outdir, grid)

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
import pickle

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import data

atmosphere_interpolator = None


def load_data():
    """Load up the pickled interpolator."""
    global atmosphere_interpolator

    filename = "Kurucz_grid_interpolator.pickle"
    atmosphere_interpolator = pickle.load(pkg_resources.open_binary(data,
                                                                    filename))


def atmos(teff, logg, feh):
    """Calculate the interpolated atmosphere for a set of stellar parameters.

    Parameters
    ----------
    teff : float
        Effective temperature of the star
    logg : float
        Surface gravity of the star
    feh : float
        [Fe/H] of the star

    Returns
    -------
    model_atmosphere : list of [RHOX,T,P,XNE]
        Provides the interpolated model atmosphere

    """
    global atmosphere_interpolator

    if atmosphere_interpolator is None:
        load_data()

    points = [[teff, logg, feh, x] for x in np.arange(72)]

    return atmosphere_interpolator(points)


def print_output(output, outfile, teff, logg, feh, vt):
    """Print the interpolated atmosphere to a file.

    Parameters
    ----------
    output : float
        List of [RHOX,T,P,XNE] describing model atmosphere
    outfile : string
        Name of the output file
    teff : float
        Effective temperature of the star
    logg : float
        Surface gravity of the star
    feh : float
        [Fe/H] of the star
    vt : float
        Velocity of microturbulence of the star

    """
    vt *= 1e5

    with open(outfile, 'w') as f:
        f.write("KURUCZ\n")
        f.write("#OVER72: T= {:.0f},[g]={:.2f},[Fe/H]={:.2f},vt={:.2e}\n".format(teff,logg,feh,vt))
        f.write("NTAU            72\n")

        for i, line in enumerate(output):
            f.write(" {:.8e}   {:.1f} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}\n".format(line[0],
                                                                                   line[1],
                                                                                   line[2],
                                                                                   line[3],
                                                                                   line[4],
                                                                                   line[5],
                                                                                   line[6]))

        f.write("     {:.2e}\n".format(vt))
        f.write("NATOMS           1  {:.2f}\n".format(feh)) #these should be editable to meet users needs
        f.write("      3.00     3.30\n")
        f.write("NMOL            19\n")
        f.write("     607.0     108.0     106.0     107.0\n")
        f.write("     606.0     608.0     112.0     707.0\n")
        f.write("     708.0     808.0      12.1   60808.0\n")
        f.write("   10108.0     101.0       6.1       7.1\n")
        f.write("       8.1     822.0      22.1\n")



def create_kurucz_array():

    data_dir = "/usr/local/mspawn/"

    files = ['AM01K2.DAT', 'AM02K2.DAT', 'AM03K2.DAT', 'AM05K2.DAT', 'AM10K2.DAT',
             'AM15K2.DAT', 'AM20K2.DAT', 'AM25K2.DAT', 'AM30K2.DAT', 'AM35K2.DAT',
             'AM40K2.DAT', 'AP00K2.DAT', 'AP01K2.DAT', 'AP02K2.DAT', 'AP03K2.DAT',
             'AP05K2.DAT']

    files = [data_dir+f for f in files]

    feh_list = [-0.1, -0.2, -0.3, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, 0.0, 0.1, 0.2, 0.3, 0.5]

    dtype = [('rho_X','f8'), ('T','f8'), ('P','f8'), ('rho_e','f8'),
             ('kappa','f8'), ('rad_acc','f8'), ('V_macro','f8')]

    full_dtype = [('rho_X','f8'), ('T','f8'), ('P','f8'), ('rho_e','f8'),
                  ('kappa','f8'), ('rad_acc','f8'), ('V_macro','f8'),
                  ('teff','f8'), ('logg','f8'), ('feh','f8'), ('step','f8')]

    teff_set = []
    logg_set = []
    feh_set = []
    data_set = []


    for j, file in enumerate(files):
        with open(file, "r") as f:

            while True:

                line = f.readline()

                # T_eff
                while "TEFF" not in line:
                    line = f.readline()
                    if line.strip() == '': break
                if line.strip() == '': break
                array = line.split()
                teff = array[1]
                logg = array[3]
                feh = feh_list[j]
                # teff_set.append(teff)
                # logg_set.append(logg)
                # feh_set.append(feh)

                # Parse through the data
                while "READ" not in line:
                    line = f.readline()
                    if line.strip() == '': break
                if line.strip() == '': break

                data_array = np.array([], dtype=dtype)

                found = False
                while True:
                    line = f.readline()
                    data = line.split()
                    if len(data) != 9 and found == True:
                        teff_set.append(teff)
                        logg_set.append(logg)
                        feh_set.append(feh)
                        data_set.append(data_array)
                        break
                    elif len(data) != 9:
                        break
                    else:
                        found = True
                    data = np.core.records.fromarrays(data,
                                                      names='rho_X, T, P, rho_e, kappa, rad_acc, V_macro, col8, col9',
                                                      formats = 'f8, f8, f8, f8, f8, f8, f8, f8, f8')
                    data = data[['rho_X', 'T', 'P', 'rho_e', 'kappa', 'rad_acc', 'V_macro']]
                    data_array = np.append(data_array, data)

                while True:
                    line = f.readline()
                    data = line.split()
                    if len(data) != 7: break
                    # data = np.core.records.fromarrays(data,
                    #                                   names='rho_X, T, P, rho_e, kappa, rad_acc, V_macro',
                    #                                   formats = 'f8, f8, f8, f8, f8, f8, f8')
                    # data_array = np.append(data_array, data)

                if line.strip() == '': break

                # data_set.append(data_array)


    full_data = np.array([], dtype=full_dtype)

    for i in range(len(teff_set)):

        new_array = np.zeros(len(data_set[i]), dtype=full_dtype)
        for col in data_set[i].dtype.names:
            new_array[col] = data_set[i][col]
        new_array['teff'] = teff_set[i]
        new_array['logg'] = logg_set[i]
        new_array['feh'] = feh_set[i]
        new_array['step'] = np.arange(len(new_array))

        full_data = np.append(full_data, new_array)

    return full_data


def create_kurucz_pickle(filename_npy):
    """Create pickle file for Kurucz data interpolation."""
    # Load Kurucz structured numpy array

    filename_kurucz_array = "kurucz_data.npy"
    kurucz_array = np.load(pkg_resources.open_binary(data,
                                                     filename_kurucz_array))

    # Convert to 2D numpy array
    atmosphere_data = kurucz_array[['rho_X',
                                    'T',
                                    'P',
                                    'rho_e',
                                    'kappa',
                                    'rad_acc',
                                    'V_macro']].copy().view((float,
                                                             len(kurucz_array.dtype.names)))

    # Create list of unique teff, logg, and feh values
    teff_set = np.unique(kurucz_array['teff'])
    logg_set = np.unique(kurucz_array['logg'])
    feh_set = np.unique(kurucz_array['feh'])

    # Create empty array of interpolation values
    interpolation_array = np.zeros((len(teff_set),
                                    len(logg_set),
                                    len(feh_set),
                                    72,
                                    7))

    # Cycle through teff, logg, and feh values, setting array for interpolation
    for i, teff in enumerate(teff_set):
        idx_teff = np.where(kurucz_array['teff'] == teff)[0]
        for j, logg in enumerate(logg_set):
            idx_logg = np.where(kurucz_array['logg'] == logg)[0]
            for k, feh in enumerate(feh_set):
                idx_feh = np.where(kurucz_array['feh'] == feh)[0]

                idx = np.intersect1d(idx_teff, idx_logg)
                idx = np.intersect1d(idx, idx_feh)

                # If no data, move on, default will be zeros
                if len(idx) == 0:
                    continue

                interpolation_array[i, j, k] = atmosphere_data[idx][:, 0:7]

    # Create Kurucz interpolator object
    kurucz_interpolator = RegularGridInterpolator((teff_set,
                                                   logg_set,
                                                   feh_set,
                                                   np.arange(72)),
                                                  interpolation_array)

    return kurucz_interpolator

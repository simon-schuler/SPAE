import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import LinearNDInterpolator

import SPAE

atmosphere_data = None
teff_set = None
logg_set = None
feh_set = None


def load_data():
    global atmosphere_data
    global teff_set
    global logg_set
    global feh_set

    path = Path(SPAE.__file__).parent / '../data/kurucz_data.npy'
    atmosphere_data = np.load(path)

    new_data = atmosphere_data.view(np.float64).reshape(atmosphere_data.shape + (-1,))

    points = new_data[:,7:10]
    points = np.unique(points, axis=0)

    teff_set = np.unique(points[:,0])
    logg_set = np.unique(points[:,1])
    feh_set = np.unique(points[:,2])



def atmos(teff, logg, feh):
    global atmosphere_data
    global teff_set
    global logg_set
    global feh_set

    if atmosphere_data is None: load_data()


    # Find closest points
    idx = np.argsort(np.abs(teff - teff_set))
    teff_tmp = np.sort(teff_set[idx[0:2]])

    idx = np.argsort(np.abs(logg - logg_set))
    logg_tmp = np.sort(logg_set[idx[0:2]])

    idx = np.argsort(np.abs(feh - feh_set))
    feh_tmp = np.sort(feh_set[idx[0:2]])

    # Create a small data set
    small_data = np.array([], dtype=atmosphere_data.dtype)

    for teff_test in teff_tmp:
        for logg_test in logg_tmp:
            for feh_test in feh_tmp:
                idx = np.where(atmosphere_data['teff'] == teff_test)[0]
                idx = np.intersect1d(idx, np.where(atmosphere_data['logg'] == logg_test)[0])
                idx = np.intersect1d(idx, np.where(atmosphere_data['feh'] == feh_test)[0])

                small_data = np.append(small_data, atmosphere_data[idx])

    # Run Interpolation
    new_data = small_data.view(np.float64).reshape(small_data.shape + (-1,))

    points = new_data[:,7:11]
    values = new_data[:,0:7]

    output_data = find_fit(teff, logg, feh, points, values)

    return output_data


def find_fit(teff, logg, feh, points, values):

    linInter = LinearNDInterpolator(points, values, fill_value=0)

    output_data = np.zeros((72, 7))

    for i in range(72):
        val = linInter((teff, logg, feh, i))
        output_data[i] = linInter((teff, logg, feh, i))

    return output_data


def print_output(output, outfile, teff, logg, feh, vt):

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

"""Misc I/O helpers: loading a pickled Spectrum_Data object, and replacing
placeholder values in a per-order array."""

import pickle
import numpy as np


def load_object(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def replace_w(array, replace_value, replace_with = 'med' ):
    """
        Use to replace individual values with either median or average of the order
        or with a specific value
        array format: [[order1],[order2],...] orderN = [x1,x2,x3,...]
        replace_value - value to be replaced
        repace_with - 'med', 'avg', or value to replace with
    """

    for i in range(len(array)):
        gd_0 = np.where(array[i] == replace_value)
        gd = np.where(array[i] != replace_value)

        for j in gd_0:
            if replace_with == 'med':
                array[i][j] = np.median(array[i][gd])
            elif replace_with == 'avg':
                array[i][j] = np.average(array[i][gd])
            else:
                array[i][j] = replace_with

    return None

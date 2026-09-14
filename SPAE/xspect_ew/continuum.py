"""Continuum point selection, used by Spectrum_Data.normalize() before the
george Gaussian Process continuum fit."""

import numpy as np
import matplotlib.pyplot as plt


class Continuum_scan():
    '''Selects points at the continuum
    '''
    def __init__(self, distx, depth):
        #values currently viewed for selection
        self.select_window = None
        #standard deviation of selected window
        #self.current_sig = None
        #Size of selection box in x axis
        self.distx = distx
        self.points_in_window = None
        #Input spectra
        self.data = None
        #Points selected as part of the continuum
        self.select_points = None
        #Relates to how deeply to move selection box into data
        self.depth = depth
        return None

    def load_data(self,x,y):
        self.data = np.array([x,y])
        self.select_points = np.zeros(len(x))
        self.points_in_window = len(self.data[0][np.where(self.data[0] <= self.data[0][0]+self.distx)])
        return None

    def scan(self):
        split_order_into = int(np.ceil(len(self.data[0])/self.points_in_window))
        split_order_x = np.array_split(self.data[0], split_order_into)
        split_order_y = np.array_split(self.data[1], split_order_into)
        for i in range(len(split_order_y)):
            dex = np.where((self.data[0] >= split_order_x[i][0])&(self.data[0] <= split_order_x[i][-1]))
            percent = np.percentile(split_order_y[i], self.depth)
            self.select_points[dex] = (split_order_y[i] >= percent)
        return None

    def view_selected(self):
        fig = plt.figure(figsize=(15,5))
        ax = fig.add_subplot(111)
        ax.scatter(self.data[0], self.data[1], c = '#cccccc', alpha = 0.75, s = 5)
        bool_points = (self.select_points == 1)
        ax.scatter(self.data[0][bool_points],self.data[1][bool_points], c = 'g', s = 5)
        plt.show()
        return None

    def get_selected(self):
        bool_points = (self.select_points == 1)
        return bool_points

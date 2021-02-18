import numpy as np

#Finding 1-sigma confidence levels
def sigma_func(x):
    x_sorted = np.sort(x)
    x_length = len(x_sorted)
    
    x_median = x_sorted[int(x_length/2)]
    x_lower = x_sorted[int(x_length*0.159)]
    x_higher = x_sorted[int(x_length*0.841)]
    
    return x_median, x_lower - x_median, x_higher - x_median



import numpy as np
from scipy.interpolate import interp1d


class PartitionFunction:
    """A class that provides the partition function for atomic elements."""

    def __init__(self):
        """Initialize a PartitionFunction object."""

        def get_T_data(array):
            """Get the numerical data from the partition function arrays."""
            dtype = array.dtype
            data = array[list(dtype.names[1:])]

            arr = np.array([data[field] for field in data.dtype.names])

            return arr

        # Load up data tables from Barklem & Collet (2016)
        dtype = [('element', 'U8')]
        for i in range(42):
            dtype.append(('T'+str(i+1), 'f8'))
        partition_filename = "../SPAE/data/J_A+A_588_A96_table8.dat.txt"
        partition_function_data = np.genfromtxt(partition_filename,
                                                delimiter='|',
                                                skip_header=8,
                                                skip_footer=1,
                                                dtype=dtype)

        dtype = [('N', 'i8'), ('element', 'U2'), ('IE1', 'f8'), ('IE2', 'f8'),
                 ('IE3', 'f8')]
        energy_filename = "../SPAE/data/J_A+A_588_A96_table4.dat.txt"
        energy_data = np.genfromtxt(energy_filename,
                                    delimiter='|',
                                    skip_header=7,
                                    skip_footer=1,
                                    dtype=dtype)

        # Temperature array from Barklem & Collet (2016)
        Temp_set = np.array([1.00000e-05, 1.00000e-04, 1.00000e-03,
                             1.00000e-02, 1.00000e-01, 1.50000e-01,
                             2.00000e-01, 3.00000e-01, 5.00000e-01,
                             7.00000e-01, 1.00000e+00, 1.30000e+00,
                             1.70000e+00, 2.00000e+00, 3.00000e+00,
                             5.00000e+00, 7.00000e+00, 1.00000e+01,
                             1.50000e+01, 2.00000e+01, 3.00000e+01,
                             5.00000e+01, 7.00000e+01, 1.00000e+02,
                             1.30000e+02, 1.70000e+02, 2.00000e+02,
                             2.50000e+02, 3.00000e+02, 5.00000e+02,
                             7.00000e+02, 1.00000e+03, 1.50000e+03,
                             2.00000e+03, 3.00000e+03, 4.00000e+03,
                             5.00000e+03, 6.00000e+03, 7.00000e+03,
                             8.00000e+03, 9.00000e+03, 1.00000e+04])

        # Create a list of elements and associated interpolation functions
        # (with respect to temperature)
        self.elements = []
        self.interp_functions = []
        for i in range(len(partition_function_data)):
            self.elements.append(partition_function_data[i]['element'])
            interpolation = interp1d(Temp_set,
                                     get_T_data(partition_function_data[i]))
            self.interp_functions.append(interpolation)

        # Convert to numpy array
        self.elements = np.array(self.elements)

    def __call__(self, element, temperature):
        """Evaluate the partition function of an element at a temperature."""
        # Error checking
        if element not in self.elements:
            raise ValueError("The requested element must be in the list "
                             "provided by Barklem & Collet (2016). Available "
                             "options are:" + str(self.elements))

        if hasattr(temperature, '__iter__'):
            if np.any(temperature < 1.0e-5) or np.any(temperature > 1.0e4):
                raise ValueError("At least one of the temperatures you have "
                                 "inputed exceed the temperature bounds (1e-5 "
                                 "K, 1e4 K) set by the partition function "
                                 "bounds from provided by Barklem & Collet "
                                 "(2016)")
        else:
            if temperature < 1.0e-5 or temperature > 1.0e4:
                raise ValueError("Your temperature exceeds the temperature "
                                 "bounds set by the partition function bounds "
                                 "(1e-5 K, 1e4 K) from provided by Barklem & "
                                 "Collet (2016)")

        # Find the index of the correct partition function
        idx = np.where(self.elements == element)[0][0]

        # Evaulate the partition function
        partition_function = self.interp_functions[idx](temperature)

        return partition_function

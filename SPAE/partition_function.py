import numpy as np
from scipy.interpolate import interp1d


class AtomicData:
    """A class that provides atomic data."""

    def __init__(self):
        """Initialize the atomic data class."""
        # Initialize the partition function data
        self.initialize_partition_function_data()

        # Initialize the ionization energy data
        self.initialize_energy_data()

    def initialize_energy_data(self):
        """Initialize the energy data."""
        # Load up data table from Barklem & Collet (2016)
        dtype = [('N', 'i8'), ('element', 'U2'), ('IE1', 'f8'), ('IE2', 'f8'),
                 ('IE3', 'f8')]
        energy_filename = "../SPAE/data/J_A+A_588_A96_table4.dat.txt"
        self.energy_data = np.genfromtxt(energy_filename,
                                         delimiter='|',
                                         skip_header=7,
                                         skip_footer=1,
                                         dtype=dtype)

        # Cycle through energy_data elements to strip white spaces
        for i in range(len(self.energy_data)):
            self.energy_data[i]['element'] = \
                self.energy_data[i]['element'].strip()

    def initialize_partition_function_data(self):
        """Initialize the partition function data for atoms."""

        def get_T_data(array):
            """Get the numerical data from the partition function arrays."""
            dtype = array.dtype
            data = array[list(dtype.names[1:])]

            arr = np.array([data[field] for field in data.dtype.names])

            return arr

        # Load up data table from Barklem & Collet (2016)
        dtype = [('element', 'U8')]
        for i in range(42):
            dtype.append(('T'+str(i+1), 'f8'))
        partition_filename = "../SPAE/data/J_A+A_588_A96_table8.dat.txt"
        partition_function_data = np.genfromtxt(partition_filename,
                                                delimiter='|',
                                                skip_header=8,
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
        self.partition_elements = []
        self.partition_interp_functions = []
        for i in range(len(partition_function_data)):
            self.partition_elements.append(
                        partition_function_data[i]['element'])
            interpolation = interp1d(Temp_set,
                                     get_T_data(partition_function_data[i]))
            self.partition_interp_functions.append(interpolation)

        # Convert to numpy array
        self.partition_elements = np.array(self.partition_elements)

    def partition_function(self, element, temperature):
        """Evaluate the partition function of an element at a temperature."""
        # Error checking
        if element not in self.partition_elements:
            raise ValueError("The requested element must be in the list "
                             "provided by Barklem & Collet (2016). Available "
                             "options are:" + str(self.partition_elements))

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
        idx = np.where(self.partition_elements == element)[0][0]

        # Evaulate the partition function
        partition_function = self.partition_interp_functions[idx](temperature)

        return partition_function

    def energy_states(self, element, ionization_state):
        """Provide the energy states for different ionization levels."""
        # Error checking
        if element not in self.energy_data['element']:
            raise ValueError("The requested element must be in the list "
                             "provided by Barklem & Collet (2016). Available "
                             "options are:"
                             + str(self.energy_data['element']))

        if ionization_state not in [1, 2, 3]:
            raise ValueError("We can only provide ionization state data for "
                             "either the first, second, or third electron.")

        # Find the correct data entry for the input element
        idx = np.where(self.energy_data['element'] == element)[0][0]

        # Create the array
        IE_list = ['IE1', 'IE2', 'IE3']

        # Find and return the requested ionization energy
        ionization_energy = self.energy_data[idx][IE_list[ionization_state-1]]

        return ionization_energy

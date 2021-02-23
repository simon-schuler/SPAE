#Writing MOOG parameter file for the parameter, abundance, and error calculations.
#The parameter file only needs to be written once, at beginning of the routine, because the output 
#files are overwritten with each itereation of the routine, only minimal output data are needed.
#
#The user can choose to have the parameter file written to screen by choosing verbose=True
#The user can choose to have more detailed MOOG output by chooseing the appropriate values for the
#MOOG input parameters.
import numpy as np

def param_file(linelist,atmosphere=0,molecules=1,lines=0,flux=0,damp=0,plot=0,units=0,verbose=False):
    if verbose:
        print('abfind')
        print('terminal        \'x11\'')
        print('standard_out    \'moog_out.1\'')
        print('summary_out     \'moog_out.2\'')
        print('model_in        \'star.mod\'')
        print('lines_in        \'' + linelist + '\'')
        print('atmosphere     ' + str(atmosphere))
        print('molecules      ' + str(molecules))
        print('lines          ' + str(lines))
        print('flux/int       ' + str(flux))
        print('damping        ' + str(damp))
        print('plot           ' + str(plot))
        print('units          ' + str(units))

    with open('batch.par', 'wt') as file:
        file.write('abfind' + '\n')
        file.write('terminal        \'x11\'' + '\n')
        file.write('standard_out    \'moog_out.1\'' + '\n')
        file.write('summary_out     \'moog_out.2\'' + '\n')
        file.write('model_in        \'star.mod\'' + '\n')
        file.write('lines_in        \'' + linelist + '\'' + '\n')
        file.write('atmosphere     ' + str(atmosphere) + '\n')
        file.write('molecules      ' + str(molecules) + '\n')
        file.write('lines          ' + str(lines) + '\n')
        file.write('flux/int       ' + str(flux) + '\n')
        file.write('damping        ' + str(damp) + '\n')
        file.write('plot           ' + str(plot) + '\n')
        file.write('units          ' + str(units) + '\n')


#Function for creating the solar and stellar linelists
def linelist_create(star_in, sun_in, direc_path):
    
    with open(direc_path + '/linelist_star.txt', 'w') as out_star: 
        with open(direc_path + '/linelist_sun.txt', 'w') as out_sun: 
            with open(star_in) as file_star: 
                with open(sun_in) as file_sun: 
                    
                    line_star = file_star.readline()
                    out_star.write(line_star) #accounts for comment line in linelist files
                    
                    line_sun = file_sun.readline()
                    out_sun.write(line_sun) #accounts for comment line in linelist files
                    
                    line = file_star.readlines()
                    line_s = file_sun.readlines()
                    
                    for line_star in line:
                        line_star_split = line_star.split()
                        #if len(line_star_split) < 2: continue
                        for line_sun in line_s:
                            line_sun_split = line_sun.split()
                            #if len(line_sun_split) < 2: continue
                            if line_star_split[0] == line_sun_split[0] and line_star_split[1] == line_sun_split[1]:
                                out_star.write(line_star)
                                out_sun.write(line_sun)
                                continue


#Reads Moog output files, parsing elements and colums
def read_file(filename):
    count = 0
    elements = ['Fe I ', 'Fe II ', 'C  I ', 'N  I ', 'O  I ', 'S  I', 'K  I ', 'Na I ', 'Mg I ', 'Al I ', 'Si I ', 'Ca I ', 'Sc II ', 'Ti I ', 'Ti II ', 'V  ', 'Cr I ',
                'Mn I ', 'Co I ', 'Ni I ', 'Cu I ', 'Zn I ', 'Ba II ']

    dtype = [('wavelength', 'f8'),
             ('ID', 'f8'),
             ('EP', 'f8'),
             ('logGF', 'f8'),
             ('EWin', 'f8'),
             ('logRWin', 'f8'),
             ('abund', 'f8'),
             ('delavg', 'f8')]

    abundances = []
    el_found = []

    with open(filename) as file:
        while True:
            count += 1

            # Get next line from file
            line = file.readline()

            # if line is empty end of file is reached
            if not line: break

            for j, el in enumerate(elements):
                species = 'Abundance Results for Species ' + el
                if species in line:
                    new_arr = []
                    el_found.append(el)
                    line = file.readline().split()
                    line = file.readline().split()
                    
                    while len(line) == 8:
                        new_arr.append(line)
                        line = file.readline().rstrip().split()
                    
                    new_arr = np.array(new_arr)
                    new_arr = np.core.records.fromarrays(new_arr.T,dtype=dtype)
                    


                    abundances.append(new_arr)

                    
    return el_found, abundances



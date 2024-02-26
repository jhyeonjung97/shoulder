# Ref: https://vasppy.readthedocs.io/en/latest/_modules/vasppy/doscar.html

import time
import datetime
import argparse
import numpy as np
import pandas as pd
from scipy.integrate import simpson

parser = argparse.ArgumentParser(description='Command-line options example')

parser.add_argument('-a', '--atoms', type=str, default='', help='atoms (e.g. "14,15", "14-16")')
parser.add_argument('-e', '--energy', type=str, default='', help='energy range (e.g., "-10,5")')
parser.add_argument('-o', '--orbital', type=str, default='d', help='orbital (e.g. "s", "p", "d", "f"')
parser.add_argument('-m', '--subset', type=str, default='', help='orbital subset # 1 for d_xy, 2 for d_yz, 3 for d_z2-r2, 4 for d_xz, 5 for d_x2-y2)')

# args, remaining_args = parser.parse_known_args()
args = parser.parse_args()
        
# Process arguments parsed by argparse
orb = args.orbital

if args.energy:
    emin, emax = map(int, args.energy.split(','))
else:
    emin, emax = None, None

# Check if input_str contains a dash, indicating a range
if '-' in args.atoms:
    start, end = args.atoms.split('-')
    atoms = list(range(int(start), int(end)+1))
    print(atoms)
elif args.atoms:
    atoms = list(map(int, args.atoms.split(',')))
    print(atoms)

if args.subset:
    m = []
    subset_numbers = list(map(int, args.subset.split(',')))
    if orb == 'p':
        subset_dict = {
            1: 'x',
            2: 'y',
            3: 'z'
        }
        m = [subset_dict[number] for number in subset_numbers if number in subset_dict]
    elif orb == 'd':
        subset_dict = {
            1: 'xy',
            2: 'yz',
            3: 'z2-r2',
            4: 'xz',
            5: 'x2-y2'
        }
        m = [subset_dict[number] for number in subset_numbers if number in subset_dict]
    elif orb == 'f':
        subset_dict = {
            1: 'y(3x2-y2)', 
            2: 'xyz', 
            3: 'yz2', 
            4: 'z3', 
            5: 'xz2', 
            6: 'z(x2-y2)', 
            7: 'x(x2-3y2)'
        }
        m = [subset_dict[number] for number in subset_numbers if number in subset_dict]
    print(m)

def pdos_column_names(lmax, ispin):
    if lmax == 1:
        names = ['s', 'p_y', 'p_z', 'p_x']
    elif lmax == 2:
        names = ['s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2']
    elif lmax == 3:
        names = ['s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2',
                  'f_y(3x2-y2)', 'f_xyz', 'f_yz2', 'f_z3', 'f_xz2', 'f_z(x2-y2)', 'f_x(x2-3y2)']
    else:
        raise ValueError('lmax value not supported')
    if ispin == 2:
        all_names = []
        for n in names:
            all_names.extend(['{}_up'.format(n), '{}_down'.format(n)])
    else:
        all_names = names
    all_names.insert(0, 'energy')
    return all_names

class Doscar:
    '''
    Contains all the data in a VASP DOSCAR file, and methods for manipulating this.
    '''

    number_of_header_lines = 6

    def __init__(self, filename, ispin=2, lmax=2, lorbit=11, spin_orbit_coupling=False, read_pdos=True, species=None):
        '''
        Create a Doscar object from a VASP DOSCAR file.
        Args:
            filename (str): Filename of the VASP DOSCAR file to read.
            ispin (optional:int): ISPIN flag. 
                Set to 1 for non-spin-polarised or 2 for spin-polarised calculations.
                Default = 2.
            lmax (optional:int): Maximum l angular momentum. (d=2, f=3). Default = 2.
            lorbit (optional:int): The VASP LORBIT flag. (Default=11).
            spin_orbit_coupling (optional:bool): Spin-orbit coupling (Default=False).
            read_pdos (optional:bool): Set to True to read the atom-projected density of states (Default=True).
            species (optional:list(str)): List of atomic species strings, e.g. ['Fe', 'Fe', 'O', 'O', 'O'].
                Default=None.
        '''
        self.filename = filename
        self.ispin = ispin
        self.lmax = lmax
        self.spin_orbit_coupling = spin_orbit_coupling
        if self.spin_orbit_coupling:
            raise NotImplementedError('Spin-orbit coupling is not yet implemented')
        self.lorbit = lorbit
        self.pdos = None
        self.species = species
        self.read_header()
        self.read_total_dos()
        if read_pdos:
            try:
                self.read_projected_dos()
            except:
                raise
        # if species is set, should check that this is consistent with the number of entries in the
        # projected_dos dataset
        
    @property
    def number_of_channels(self):
        if self.lorbit == 11:
            return {1:4, 2: 9, 3: 16}[self.lmax]
        raise NotImplementedError

    def read_header(self):
        self.header = []
        with open(self.filename, 'r') as file_in:
            for i in range(Doscar.number_of_header_lines):
                self.header.append(file_in.readline())
        self.process_header()

    def process_header(self):
        self.number_of_atoms = int(self.header[0].split()[0])
        self.number_of_data_points = int(self.header[5].split()[2])
        self.efermi = float(self.header[5].split()[3])
        
    def read_total_dos(self): # assumes spin_polarised
        start_to_read = Doscar.number_of_header_lines
        df = pd.read_csv(self.filename, 
                         skiprows=start_to_read, 
                         nrows=self.number_of_data_points,
                         delim_whitespace=True, 
                         names=['energy', 'up', 'down', 'int_up', 'int_down'],
                         index_col=False)
        self.energy = df.energy.values
        df.drop('energy', axis=1)
        self.tdos = df
        
    def read_atomic_dos_as_df(self, atom_number): # currently assume spin-polarised, no-SO-coupling, no f-states
        assert atom_number > 0 & atom_number <= self.number_of_atoms
        start_to_read = Doscar.number_of_header_lines + atom_number * (self.number_of_data_points + 1)
        df = pd.read_csv(self.filename,
                         skiprows=start_to_read,
                         nrows=self.number_of_data_points,
                         delim_whitespace=True,
                         names=pdos_column_names(lmax=self.lmax, ispin=self.ispin),
                         index_col=False)
        return df.drop('energy', axis=1)
    
    def read_projected_dos(self):
        """
        Read the projected density of states data into """
        pdos_list = []
        for i in range(self.number_of_atoms):
            df = self.read_atomic_dos_as_df(i+1)
            pdos_list.append(df)
        self.pdos = np.vstack([np.array(df) for df in pdos_list]).reshape(
            self.number_of_atoms, self.number_of_data_points, self.number_of_channels, self.ispin)
        
    def pdos_select(self, atoms=None, spin=None, l=None, m=None):
        """
        Returns a subset of the projected density of states array.
        """
        valid_m_values = {'s': [],
                          'p': ['x', 'y', 'z'],
                          'd': ['xy', 'yz', 'z2-r2', 'xz', 'x2-y2'],
                          'd1': ['xy'],
                          'd2': ['yz'],
                          'd3': ['z2-r2'],
                          'd4': ['xz'],
                          'd5': ['x2-y2'],
                          'f': ['y(3x2-y2)', 'xyz', 'yz2', 'z3', 'xz2', 'z(x2-y2)', 'x(x2-3y2)']}
        if not atoms:
            atom_idx = list(range(self.number_of_atoms))
        else:
            atom_idx = atoms
        to_return = self.pdos[atom_idx, :, :, :]
        if not spin:
            spin_idx = list(range(self.ispin))
        elif spin == 'up':
            spin_idx = [0]
        elif spin == 'down':
            spin_idx = [1]
        # elif spin == 'both':
        #     spin_idx = [0,1]
        else:
            raise ValueError 
        to_return = to_return[:, :, :, spin_idx]
        if not l:
            channel_idx = list(range(self.number_of_channels))
        elif l == 's':
            channel_idx = [0]
        elif l == 'p':
            if not m:
                channel_idx = [1, 2, 3]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['p']) if v in m]
        elif l == 'd':
            if not m:
                channel_idx = [4, 5, 6, 7, 8]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['d']) if v in m]
        elif l == 'f':
            if not m:
                channel_idx = [9, 10, 11, 12, 13, 14, 15]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['f']) if v in m]
        else:
            raise ValueError
        to_return = to_return[:, :, channel_idx, :]
        return to_return
    
    def pdos_sum(self, atoms=None, spin=None, l=None, m=None):
        return np.sum(self.pdos_select(atoms=atoms, spin=spin, l=l, m=m), axis=(0,2,3))

if 'f' in orb:
    lmax = 3
elif 'd' in orb:
    lmax = 2
elif 'p' in orb:
    lmax = 1
else:
    print('lmax value not supported')
with open('OUTCAR', 'r') as file:
    for line in file:
        if 'ISPIN' in line and '1' in line:
            ispin = 1
        elif 'ISPIN' in line and '2' in line:
            ispin = 2

# calculation of d-band center
# Open doscar
dosfile = 'DOSCAR'
doscar  = Doscar(dosfile, ispin=ispin, lmax=lmax, lorbit=11)  # calculation setting 

# Set atoms for integration
if ispin == 1:
    non = doscar.pdos_sum(atoms, spin='up', l=orb)
elif ispin == 2:
    up = doscar.pdos_sum(atoms, spin='up', l=orb)
    down = doscar.pdos_sum(atoms, spin='down', l=orb)
else:
    print('ispin value not supported')    

# Set intergrating range 
efermi = doscar.efermi #- doscar.efermi 
energies = doscar.energy - doscar.efermi
if emin == None:
    emin = energies[0]
if emax == None:
    emax = energies[-1]
erange = (emin, emax)
emask = (energies <= erange[-1])

# Calculating center of the orbital specified above in line 184
x = energies[emask]
if ispin == 1:
    y = non[emask]
    dbc = simpson(y=y*x, x=x) / simpson(y=y, x=x)
    print('dbc: {:.4f} (eV)'.format(dbc))
elif ispin == 2:
    y1 = up[emask]
    y2 = down[emask]
    dbc_up   = simpson(y=y1*x, x=x) / simpson(y=y1, x=x)
    dbc_down = simpson(y=y2*x, x=x) / simpson(y=y2,x=x)
    dbc = simpson(y=(y1+y2)*x, x=x) / simpson(y=(y1+y2), x=x)
    print('dbc_up  : {:.4f} (eV)'.format(dbc_up))
    print('dbc_down: {:.4f} (eV)'.format(dbc_down))
    print('dbc     : {:.4f} (eV)'.format(dbc))

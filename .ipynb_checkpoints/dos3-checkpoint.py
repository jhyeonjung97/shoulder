# Ref: https://vasppy.readthedocs.io/en/latest/_modules/vasppy/doscar.html

import time
import datetime
import argparse
import numpy as np
import pandas as pd
from scipy.integrate import simpson

parser = argparse.ArgumentParser(description='Command-line options example')

parser.add_argument('-a', '--atoms', type=str, default='', help='atoms (e.g. "14,15", "14-16")')
parser.add_argument('-e', '--energy', type=str, default='', help='energy range (e.g., -e=-10,5)')
parser.add_argument('-o', '--orbital', type=str, default='d', help='orbital (e.g. "s", "p", "d", "f"')
parser.add_argument('-m', '--subset', type=str, default='', help='orbital subset # 1 for d_xy, 2 for d_yz, 3 for d_z2-r2, 4 for d_xz, 5 for d_x2-y2)')

# args, remaining_args = parser.parse_known_args()
args = parser.parse_args()
        
# Process arguments parsed by argparse
orb = args.orbital

if args.energy:
    emin, emax = map(float, args.energy.split(','))
else:
    emin, emax = None, None

# Check if input_str contains a dash, indicating a range
if '-' in args.atoms:
    a_start, a_end = args.atoms.split('-')
    atoms = list(range(int(a_start), int(a_end)+1))
elif args.atoms:
    atoms = list(map(int, args.atoms.split(',')))

if args.subset:
    if '-' in args.subset:
        m_start, m_end = args.subset.split('-')
        subset_numbers = list(range(int(m_start), int(m_end)+1))
    else:
        subset_numbers = list(map(int, args.subset.split(',')))
    if orb == 's':
        subset_dict = {
            1: ' '
        }
        m = [subset_dict[number] for number in subset_numbers if number in subset_dict]
    elif orb == 'p':
        subset_dict = {
            1: 'y',
            2: 'z',
            3: 'x'
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
    print("\033[1m\033[34m" + ' '.join([orb + '_' + mm for mm in m]) + "\033[0m")
else:
    if orb == 's':
        print('s')
    m = None

def pdos_column_names(lmax, ispin):
    if lmax == 0:
        names = ['s']
    elif lmax == 1:
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
        
    @property
    def number_of_channels(self):
        if self.lorbit == 11:
            return {0: 1, 1: 4, 2: 9, 3: 16}[self.lmax]
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
        valid_m_values = {'s': [' '],
                          'p': ['y', 'z', 'x'],
                          'd': ['xy', 'yz', 'z2-r2', 'xz', 'x2-y2'],
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
            if not m:
                channel_idx = [0]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['s']) if v in m]
        elif l == 'p':
            if not m:
                channel_idx = [1, 2, 3]
            else:
                channel_idx = [i+1 for i, v in enumerate(valid_m_values['p']) if v in m]
        elif l == 'd':
            if not m:
                channel_idx = [4, 5, 6, 7, 8]
            else:
                channel_idx = [i+4 for i, v in enumerate(valid_m_values['d']) if v in m]
        elif l == 'f':
            if not m:
                channel_idx = [9, 10, 11, 12, 13, 14, 15]
            else:
                channel_idx = [i+9 for i, v in enumerate(valid_m_values['f']) if v in m]
        else:
            raise ValueError
        to_return = to_return[:, :, channel_idx, :]
        return to_return, len(channel_idx)
    
    def pdos_sum(self, atoms=None, spin=None, l=None, m=None):
        pdos_subset, channel_idx_length = self.pdos_select(atoms=atoms, spin=spin, l=l, m=m)
        return np.sum(pdos_subset, axis=(0,2,3)), channel_idx_length

def check_orbitals_in_potcar(potcar_path):
    has_d_or_f_orbital = False
    orbital_types = []

lmax = None
with open('POTCAR', 'r') as file:
    for line in file:
        if "VRHFIN" in line:
            if 'f' in line:
                lmax = 3
            elif 'd' in line and lmax != 3:
                lmax = 2
            elif 'p' in line and lmax != 2:
                lmax = 1
            elif 's' in line and lmax == None:
                lmax = 0
if lmax == None:
    print('check lmax value..')

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
    non, o_num = doscar.pdos_sum(atoms, spin='up', l=orb, m=m)
elif ispin == 2:
    up, o_num_up = doscar.pdos_sum(atoms, spin='up', l=orb, m=m)
    down, o_num_down = doscar.pdos_sum(atoms, spin='down', l=orb, m=m)
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
emask_occ = (energies <= 0)
emask_unocc = (energies > 0)

# Calculating center of the orbital specified above in line 184
if ispin == 1:
    x = energies[emask]
    y = non[emask]
    dbc = simpson(y=y*x, x=x) / simpson(y=y, x=x)
    print('  dbc: {:.4f} (eV)\n'.format(dbc))
    x = energies[emask_occ]
    y = non[emask_occ]
    occ = simpson(y=y, x=x)
    x = energies[emask_unocc]
    y = non[emask_unocc]
    unocc = simpson(y=y, x=x)
    e_num = occ / unocc * o_num * 2
    print('  occ: {:.4f} (eV)\n'.format(occ),
          'unocc: {:.4f} (eV)\n'.format(unocc),
          'e_num: {:.4f} (eV)\n'.format(e_num))
elif ispin == 2:
    x = energies[emask]
    y1 = up[emask]
    y2 = down[emask]
    dbc_up   = simpson(y=y1*x, x=x) / simpson(y=y1, x=x)
    dbc_down = simpson(y=y2*x, x=x) / simpson(y=y2, x=x)
    dbc = simpson(y=(y1+y2)*x, x=x) / simpson(y=(y1+y2), x=x)
    # print('   dbc_up  : {:.4f} (eV)\n'.format(dbc_up),
    #       '  dbc_down: {:.4f} (eV)\n'.format(dbc_down),
    #       '  dbc     : {:.4f} (eV)\n'.format(dbc))
    total1 = simpson(y=y1, x=x)
    total2 = simpson(y=y2, x=x)
    total = total1+total2
    x = energies[emask_occ]
    y1 = up[emask_occ]
    y2 = down[emask_occ]
    occ1 = simpson(y=y1, x=x)
    occ2 = simpson(y=y2, x=x)
    occ = occ1+occ2
    dbc_occ_up = simpson(y=y1*x, x=x) / simpson(y=y1, x=x)
    dbc_occ_down = simpson(y=y2*x, x=x) / simpson(y=y2, x=x)
    x = energies[emask_unocc]
    y1 = up[emask_unocc]
    y2 = down[emask_unocc]
    unocc1 = simpson(y=y1, x=x)
    unocc2 = simpson(y=y2, x=x)
    unocc = unocc1+unocc2
    e_num1 = occ1 / total1 * o_num_up * 2
    e_num2 = occ2 / total2 * o_num_down * 2
    e_num = occ / total * o_num_up * 2
    dbc_unocc_up = simpson(y=y1*x, x=x) / simpson(y=y1, x=x)
    dbc_unocc_down = simpson(y=y2*x, x=x) / simpson(y=y2, x=x)
    print('   occ:\t{:.4f}\t{:.4f} (eV)\n'.format(dbc_occ_up, dbc_occ_down),
          '  occ:\t{:.4f}\t{:.4f} (e-)\n'.format(occ1, occ2),
          'unocc:\t{:.4f}\t{:.4f} (eV)\n'.format(dbc_unocc_up, dbc_unocc_down),
          'unocc:\t{:.4f}\t{:.4f} (e-)\n'.format(unocc1, unocc2),
          'total: {:.4f}\t{:.4f}(e-)\n'.format(total1, total2),
          'e_num: {:.4f}\t{:.4f}(e-)\n'.format(e_num1, e_num2))
    # print('   occ_down: {:.4f} (eV)\n'.format(dbc_occ_down),
    #       '  occ_down: {:.4f} (e-)\n'.format(occ2),
    #       'unocc_down: {:.4f} (eV)\n'.format(dbc_unocc_down),
    #       'unocc_down: {:.4f} (e-)\n'.format(unocc2),
    #       'total_down: {:.4f} (e-)\n'.format(total2),
    #       'e_num_down: {:.4f} (e-)\n'.format(e_num2))
    # print('  dbc  : {:.4f} eV'.format(dbc))
    # if not 1.9 * o_num_up * 2 < total < 2.1 * o_num_up * 2:
    #     print('  #occ : {:.4f}\n'.format(occ),
    #           '#unocc : {:.4f}\n'.format(unocc),
    #           '#total : {:.4f}'.format(total))
    # print('  occ  : {:.4f} e-'.format(e_num))
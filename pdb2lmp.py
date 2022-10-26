from  ase.io import read, write
import ase.io.vasp

atoms = read('mixture.xyz')
#pos   = atoms.get_positions()
#print(pos) 
#ase.io.vasp.write_vasp('POSCAR',atoms,label='KShim',direct=False,sort=False)
write('mixture.vasp', atoms, format='vasp') #, units='real', atom_style='charge', cube=60.)

#https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.write

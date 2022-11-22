from ase.io import read, write
from os import listdir
import subprocess
from sys import argv
import os

directory = '/Users/hailey/figure/'
select = argv[1:]
files = ' '.join(select)
here = os.getcwd()
subprocess.call('cp '+files+' '+directory+' \n cd '+directory+'', shell=True)
x, y, z = input("rotation anple (x, y, z): ").split()
x = int(x)
y = int(y)
z = int(z)
for file in listdir(directory):
    if file.endswith(".vasp") or file.endswith(".xyz"):
        filename = file.split('.')[0]
        atoms = read(file)
        atoms.center()
        write("%s/%s.png" % (directory, filename), atoms, rotation='%dx, %dy, %dz' % (x, y, z))
subprocess.call('mv '+directory+'/*.png '+here+' \n cd '+here+' \n rm '+directory+'/*', shell=True)

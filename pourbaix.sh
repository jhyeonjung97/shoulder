scp jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/1Fe_energies.tsv .
scp jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/2Co_energies.tsv .
scp jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/3Mo_energies.tsv .

python ~/bin/shoulder/pourbaix-fe.py
python ~/bin/shoulder/pourbaix-co.py
python ~/bin/shoulder/pourbaix-mo.py
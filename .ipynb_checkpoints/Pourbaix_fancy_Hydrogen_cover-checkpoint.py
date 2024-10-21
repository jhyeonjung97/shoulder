#!/usr/bin/env python

from math import pow
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.ticker import *
from scipy import *
from ase import Atoms
from ase.io import *
import subprocess

#Import plot_settings as ps
from matplotlib import rc
# with golden ration and the whole shebang ...
# settings size and font for revtex stylesheet
fig_width_pt = 1.8*246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt *= 300./72 # convert to 300 dpi
inches_per_pt = 1.0/72.27               # Convert pt to inches
#inches_per_pt = 1.0/300               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height]

font_size = 10
tick_font_size = 10
xlabel_pad = 8
ylabel_pad = 18

matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['font.serif'] = 'Arial'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.labelsize'] = 2*font_size
matplotlib.rcParams['legend.fontsize'] = font_size
matplotlib.rcParams['xtick.labelsize'] = tick_font_size
matplotlib.rcParams['ytick.labelsize'] = tick_font_size
matplotlib.rcParams['mathtext.default'] = 'regular'


matplotlib.rcParams['lines.linewidth'] = 1.
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

Umin=-0.5
Umax=2.5

#atoms=read('nanoribbon/emptyempty/moments.traj')
#for i in atoms:
#        if (i.symbol!='H' and i.symbol!='O'):
#                metal=i.symbol

#Values for X,Y
pH = arange(0,14,0.10)
U = arange(Umin,Umax,0.05)
Umax2=Umax+0.06*14
U2= arange(Umin,Umax2,0.05)

ax.axis([0, 14, Umin,Umax])
ax.set_xlabel(r'pH')
ax.set_ylabel(r'U/V')

#plt.xlabel('$pH$',family = 'Helvetica', size='40')
#plt.ylabel('$U/V$',family = 'Helvetica',size='40')

#Constants
kbt = 0.0256 
const = kbt*log(10)

kjmol = 96.485

extraticks=[1.23]
plt.yticks(list(plt.yticks()[0]) + extraticks)


def frhe(x):
    return 1.78-x*const
def fh2oequi(x):
    return 1.23-x*const

#--------------------------------------------
# Total  Energies of the different coverages |
#--------------------------------------------


surfs = [
[-962.55524065,0,0,0],# done
#[-342.310,0,-1,0],# done
[-973.60202128,2,0,0],# done 
#[-350.84760688,0,1,0],# done
#[-356.52219009,0,0,1],# done
[-979.55138349,3,0,0],#  done
#[-366.92897209,0,0,1],# done
#[-375.58742672,0,0,2],# done
[-984.74327528,4,0,0],#  done
]

#[clean surf, #Hs, #Os, #OHs] 
print (surfs)
#------------------------
# Calculation of the DG |
#------------------------

##400 eV ,O_s
#h2=  -6.759300 
#h2o= -14.01977

#500 eV, O regular
#h2=-6.77014123
#h2o=-14.21744725
h2=  -6.77149190
h2o= -14.23091949

#600 eV, O regular
#h2= -6.77412273
#h2o= -14.23797429

#Max claculated
#H2	6.770141	-	0.26810	0.09048	0.00136	0.408000	6.720721
#H2O	-14.217447	-	0.56733	0.00010	0.00186	0.558000	-14.208017

#Contributions to Gibbs Energies for gas molecules (VASP-PBE calculated by Max; T= 300K)
zpeh2o=0.560	#exp. NIST 0.558
zpeh2=0.268 	#exp. NIST 0.273
cvh2o=0.103  	#From Colin at P = 0.035 bar
cvh2=0.0905
tsh2o=0.675 	#From Colin at P = 0.035 bar
tsh2=0.408 	#at P = 1bar

#Contributions to Gibbs Energies for adsorbates (VASP-PBE calculated by Max using Michal data for NiCe; T= 300K)
#if(PBE_VASP):

zpeoh =  0.376
zpeo  =  0.064
zpeooh=  0.471
cvoh  =  0.042
cvo   =  0.034
cvooh =  0.077
tsoh  =  0.066
tso   =  0.060
tsooh =  0.134

#Gibbs Energies for the gas molecules
dgh2o=zpeh2o +cvh2o -tsh2o 
dgh2=zpeh2 +cvh2 -tsh2 
	
#Gibbs Energy corrections for adsorbates
dso=zpeo +cvo -tso -(dgh2o -dgh2)
dsoh=zpeoh +cvoh -tsoh -(dgh2o -0.5*dgh2)
dsooh=zpeooh +cvooh -tsooh -(2*dgh2o -1.5*dgh2)
dsh=dsoh-dso

#dso=zpeo -(zpeh2o-zpeh2 -tsh2o+tsh2)
#dsoh=zpeoh -(zpeh2o -0.5*zpeh2 -tsh2o+ 0.5*tsh2)
#dsooh=zpeooh-(2*zpeh2o -1.5*zpeh2 -2*tsh2o+ 1.5*tsh2)
#dsh = dsoh-dso 
print (dsoh, dso, dsh)

def addO(x,y):
 return -(h2o-h2)-2*(y+x*const)+dso
def addOH(x,y):
 return -(h2o-0.5*h2)-(y+x*const)+dsoh
def addOOH(x,y):
 return -(2*h2o-1.5*h2)-3*(y+x*const)+dsooh

def addH2O(x,y):
 return -(h2o)-(zpeh2o-tsh2o)

def addH(x,y):
 return -0.5*h2+1*(y+x*const)+dsh


print (addH(0,0), addH(0,1.4), dsh)

i=2
x=0;y=1.23
print (surfs[i][0] -surfs[0][0] +surfs[i][1]*addH(x,y) +surfs[i][2]*addO(x,y) +surfs[i][3]*addOH(x,y), surfs[i][0], surfs[0][0] , surfs[i][3], addOH(x,y))

#Function to calculate DG
def dg(i,x,y):
    return surfs[i][0] -surfs[0][0] +surfs[i][1]*addH(x,y) +surfs[i][2]*addO(x,y) +surfs[i][3]*addOH(x,y)



#find intersects

nsurfs=len(surfs)

i=0  #Ph=0
lowest_surfaces=[]
for j in U2:
    #print i,j
    values=[]
    
    for k in range(nsurfs):
            #print k,dg(k,i,j)
            values.append(dg(k,i,j))
    sorted_values=sorted(range(len(values)), key=lambda k: values[k])
    lowest_surfaces.append(sorted_values[0])
    #print j, values
    
#print lowest_surfaces

crossover=[]
uniquesurf=[]
uniquesurf.append(lowest_surfaces[0])
old_value=lowest_surfaces[0]
crossover.append(Umin)
for j in range(len(U2)):
    if(lowest_surfaces[j]!=old_value):
        uniquesurf.append(lowest_surfaces[j])
        crossover.append(U2[j])
        old_value=lowest_surfaces[j]
    
crossover.append(Umax2)    


print (crossover)
print (uniquesurf)


#color=['turquoise', 'green', 'red','blue', 'gray', 'gold', 'purple', 'pink', 'orange']
#color=['turquoise', 'green', 'red','blue', 'gray', 'gold', 'purple', 'pink', 'orange'\
#       ,'olive','yellowgreen','violet','navy','lime' ]
color=['turquoise', 'green', 'red','blue', 'gray', 'gold', 'purple', 'pink', 'darkorange'\
       ,'lime','olive','yellowgreen','violet','navy','brown','teal','deeppink',\
       'cyan','dogerblue','steelblue','darkslategrey']       
pH2 = arange(0,14,0.01)


for i in range(len(uniquesurf)):
    k=uniquesurf[i]
    foo= r"S$_{%i}$(H-%i O-%i OH-%i)" % (k, surfs[k][1], surfs[k][2],surfs[k][3])
    #fbk = {'lw':0.0, 'edgecolor':color[i]}
    fbk = {'lw':0.5, 'edgecolor':'black'}
    fill_between(pH2,crossover[i]-pH2*const, crossover[i+1]-pH2*const, facecolor=color[i],alpha=0.3,**fbk)
    plot([], [], color=color[i],alpha=0.3, linewidth=5,label=foo)

Vover=0.184
y=1.23+Vover -pH2*const
llabel='$\eta$ OER = '+ repr(Vover)+' '
plot(pH2,y,'-',color='black', lw=1, dashes=(3,1),label=llabel)

#yy=1.23-Vover_orr -pH2*const
#llabel2='$\eta$ ORR  = '+ repr(Vover_orr)+' '
#plot(pH2,yy,'-',color='green', lw=1, dashes=(3,1),label=llabel2)

plot(pH2,1.23-pH2*const,'--',color='blue',lw=1, dashes=(3,1))
ax.text(0.2,1.25,r'2H$_2$O $\leftrightarrow$ 4H$^+$ +O$_2$+4e$^-$',color='blue',rotation=-7,fontsize='x-small')
#legend() 
legend(bbox_to_anchor=(0.05, 1), loc=2, borderaxespad=0., ncol=2, fancybox=True, shadow=True, fontsize='x-small',handlelength=3)
fig.savefig('Pourbaix_2D_H_octa2.pdf', bbox_inches='tight')
fig.savefig('Pourbaix_2D_H_oct2.png', bbox_inches='tight')
#show()


fig.clf()
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1=-1.0
x2=2.5
ax.axis([-1.0, 2.5, -800, 200])


ax.set_xlabel(r'RHE (V)')
ax.set_ylabel(r'$\Delta$G (kjmol)')

xdata=[]
ydata=[]
y2data=[]


xx = np.arange(-1.0,2.5, 0.05)
#ax.plot(xx,0*xx+0,'-',lw=1, dashes=(3,1),c='black',label='OOH scaling')

for k in range(nsurfs):
    foo= r"S$_{%i}$(H: %i O: %i OH: %i)" % (k, surfs[k][1], surfs[k][2],surfs[k][3])
    ax.plot(xx, dg(k,0,xx)*kjmol, '-',lw=1, c=color[k],label=foo)


legend(bbox_to_anchor=(0.05, 1.3), loc=2, borderaxespad=0., ncol=2, fancybox=True, shadow=True, fontsize='x-small',handlelength=3)
fig.savefig('LNO_NiO_term_1D_H_octa2.pdf', bbox_inches='tight')
show()
# clear the plot
fig.clf()



 

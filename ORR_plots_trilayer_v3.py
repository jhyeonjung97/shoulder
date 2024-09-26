#! /user/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from math import pow
from pylab import *
import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.ticker import *
from scipy import *
import subprocess

#Import plot_settings as ps
from matplotlib import rc

# settings size and font for revtex stylesheet

fig_width_pt = 1.8*246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt *= 300./72 # convert to 300 dpi
inches_per_pt = 1.0/72.27               # Convert pt to inches
#inches_per_pt = 1.0/300               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height]
fig = plt.figure(figsize=fig_size,dpi=300)


font_size = 9 
tick_font_size = 8
xlabel_pad = 8
ylabel_pad = 18
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'

matplotlib.rcParams['font.size'] = 10
#matplotlib.rcParams['axes.labelsize'] = 2*font_size
matplotlib.rcParams['axes.labelsize'] = font_size
matplotlib.rcParams['legend.fontsize'] = font_size
matplotlib.rcParams['xtick.labelsize'] = tick_font_size
matplotlib.rcParams['ytick.labelsize'] = tick_font_size

#font_default='helvetica'
font_default='cmss'
#font_default='times'

def setfont(font=font_default,unicode=True):
   r"""
   Set Matplotlibs rcParams to use LaTeX for font rendering.
   Revert all changes by calling rcdefault() from matplotlib.

   Parameters:
   -----------
   font: string
       "Helvetica"
       "Times"
       "Computer Modern"

   usetex: Boolean
       Use unicode. Default: False.    
   """

   # Use TeX for all figure text!
   #plt.rc('text', usetex=True)

   font = font.lower().replace(" ","")
   if font == 'times':
       # Times
       font = {'family':'serif', 'serif':['Times']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage{mathptmx}
                   """
   elif font == 'helvetica':
       # Helvetica
       # set serif, too. Otherwise setting to times and then
       # Helvetica causes an error.
       font = {'family':'sans-serif','sans-serif':['Helvetica'],
               'serif':['cm10']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage[tx]{sfmath}
                      \usepackage{helvet}
                      \usepackage{sansmath}
                   """
   else:
       # Computer modern serif
       font = {'family':'serif', 'serif':['cm10']}
       preamble  = r"""
                      \usepackage{color}
                   """
   if font == 'cmss':
       # Computer modern sans serif
       font = {'family':'sans-serif', 'serif':['cmss']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage[tx]{sfmath}
                   """
   #if unicode:
       # Unicode for Tex
       #preamble =  r"""\usepackage[utf8]{inputenc}""" + preamble
       # inputenc should be set automatically
       #plt.rcParams['text.latex.unicode']=True
   
   #print font, preamble
   plt.rc('font',**font)
   plt.rcParams['text.latex.preamble'] = preamble



setfont(font_default,unicode=True)


ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

zoomx=0.4
zoomy=0.6
d1=4*zoomx
d2=3*zoomy
d3=5*zoomx
xcenter=0.9 #0.65
#ycenter=1.23#2.4
ycenter=3.75#2.4

x1=xcenter-d1#-0.6
x2=xcenter+d3#+0.5#2.2
y1=ycenter-d2#1#0.5
y2=ycenter+d2#+1.2#5
ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ -$\Delta$G$_{\sf O}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)')

delta = 0.01
x = np.arange(x1,x2+delta, delta)
y = np.arange(y1,y2+delta, delta)
X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)


#fit=[0.84527288, 3.38026638]
def ooh_oh_scaling(doh):
    #like ambars
    #dooh=0.5*doh  + 3.0		 #O
    #normal one
    dooh=doh + 3.2
    return dooh
def oh_ooh_scaling(dooh):
    doh=dooh=3.2
    return doh

def overpotential(doh,do):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, do-doh, dooh-do, -dooh+4.92]
    m=max(dg14)
    return m-1.23
    #return doh*do

def overpotential2(x,doh):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, x, -x+2.46, -dooh+4.92]
    m=max(dg14)
    return m-1.23
    #return doh*do

def overpotential3(x,dooh):
#    dooh=ooh_oh_scaling(doh)
    dg14=[abs(1.23-x), abs(dooh-3.69)]
    m=max(dg14)
    return m
    #return doh*do


def overpotential_orr_for_contour(doh,dooh):
    do=1.469*doh +  1.253 #from our scaling
    #do=0.21*doh +  2.06 #from our scaling, this work
    dg14=[-doh,-do+doh,-dooh+do,-4.92+dooh]
    m=max(dg14)
    return m+1.23

def overpotential_orr(doh,do,dooh):
    dg14=[-doh,-do+doh,-dooh+do,-4.92+dooh]
    m=max(dg14)
    return m+1.23

def orr_step (i):
        steps=['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']
        return str(steps[i]) 

def overpotential_orr_full(doh,do,dooh):
    #do=1.556*doh + 0.9951 #from our scaling
    dg14=[-4.92+dooh,-dooh+do,-do+doh,-doh] #add solvation soon
    m=max(dg14)
    ##print dg14.index(m)
    #print dg14,m
    print (dg14)
    return [round(m+1.23,2),round(-m,2),orr_step(dg14.index(m))]


def overpotential_label(doh,do):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, do-doh, dooh-do, dooh-4.92]
    m=max(dg14)
    for i in range(len(dg14)):
       if(m==dg14[0]):
          return r'OH lim.'
       if(m==dg14[1]):
          return r'OH-O lim.'
       if(m==dg14[2]):
          return r'O-OOH lim.'
       if(m==dg14[3]):
          return r'OOH-O$_{\sf 2}$ lim.'
    #return doh*do    


#Z=overpotential(X,Y)


#print 'ORR value ', overpotential_orr(1.23,2.46+1.23)


Z=[]
for j in y:
    tmp=[]
    for i in x:
        #tmp.append(overpotential3(i,j))
        tmp.append(overpotential_orr_for_contour(i,j))
    Z.append(tmp)

#print overpotential(0.8,2.4)

Z = np.array(Z,dtype=object)
#print Z

#im = plt.imshow(Z, origin='lower',interpolation='bilinear',
#                cmap=cm.jet_r, extent=(x1,x2,y1,y2), vmin=0, vmax=2)
import numpy as np
#import seaborn as sns
#sns.set()

from matplotlib.colors import ListedColormap


origin='lower'
levels = np.arange(0.1, 1.7, 0.1)
#levels = np.arange(0.2, 2, 0.1)
CS = plt.contourf(X, Y, Z, levels, #20, # [-1, -0.1, 0, 0.1],
                        #alpha=0.8,
                        #cmap=plt.cm.bone,
                        #cmap=plt.cm.jet_r,
                        #extend='both',
                        cmap=ListedColormap([
                        '#a50026',
                        '#d73027',
                        '#f46d43',
                        '#fdae61',
                        '#fee090',
                        '#ffffbf',
                        '#ffffe5',
                        '#ffffff',
                        '#e0f3f8',
                        '#abd9e9',
                        '#74add1',
                        '#4575b4',
                        '#313695',
                        '#d0d1e6',
                        ])  ,
                        extend='max',
                        origin=origin)

# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution,
# or leave out the levels kwarg to use all of the original levels.

'''
CS2 = plt.contour(CS, levels=CS.levels,
                        colors = 'white',
                        linewidths=0.05,
                        alpha=0.3,
                        origin=origin,
                        hold='on')
'''

#levels = np.arange(0, 2, 0.05)
#CS = plt.contourf(X,Y,Z, levels, cmap=cm.jet_r, origin='lower')
#CS = plt.contourf(X,Y,Z, levels, origin='lower')

#im = plt.imshow(Z, interpolation='bilinear', origin='lower',
#                cmap=cm.jet, extent=(x1,x2,y1,y2))


#levels2 = [2.0]
#CS2 = plt.contour(CS, levels2,
#                        colors = 'r',
#                        origin='lower',
#                        hold='on')

#CS = plt.contour(Z, levels,
#                 origin='lower',
#                 linewidths=0.5,
#                 extent=(x1,x2,y1,y2))




##Thicken the zero contour.
#zc = CS.collections[6]
#plt.setp(zc, linewidth=2)


cbar = plt.colorbar(CS, ticks=[np.arange(0.1, 1.6, step=0.1)])

#cbar.ax.set_ylabel('Overpotential [V]')
#cbar.ax.set_ylabel(r'$\eta_{\sf calc.}$')
cbar.ax.set_ylabel(r'$\eta_{\sf ORR}$ (V)')

cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')


#cbar.add_lines(CS2)


#plt.clabel(CS, levels[1::2],  # label every second level
#           inline=1,
#           fmt='%1.1f',
#           fontsize='x-small')


#plt.title('Lines with colorbar')

# We can still add a colorbar for the image, too.


# This makes the original colorbar look a bit out of place,
# so let's improve its position.

ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()


#plot(x,ooh_oh_scaling(x),'--',color='orange',lw=1, dashes=(3,1),label='$\Delta$G$_{\sf OOH}$=0.82G$_{\sf OH}$+3.18 eV')

#ax.text(x1+0.02,y2-0.3,'$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(fit[0],fit[1]),color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')

#ax.text(x1+0.02,y2-0.3,'$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(0.82,3.18),color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')


#plt.show()

offset=[0.0,0.08]


calc_systems=[] 
#energies are dGs on dEs!
#to input dEs
def make_dGs(x,i):
    dgeoh=0.339
    dgeo=0.001
    dgeooh=0.364
    #dgeoh=0.391254
    #dgeo=0.0216138
    #dgeooh=0.451538
    if(i==0):
       return x+dgeoh
    if(i==1):
       return x+dgeo
    if(i==2):
       return x+dgeooh

colors=['#a6cee3','#fdbf6f','#b2df8a','#1f78b4','#e31a1c','#fb9a99','#33a02c']
tcolor='#737373'
#tcolor='#d9d9d9'
bcolor='#252525'
#bcolor='#99000d'
symbols=['D','^','o','s','h']


# Ti
#calc_systems.append([2.714   ,  3.385 ,    5.567 r'graphene',tcolor,0.0,0.08,1.5,'8','green',symbols[4],'grey'])
#calc_systems.append([-0.582 , 0.867  , 3.965 ,1.7, r'full-hydroxy',tcolor,0.0,0.08,1.5,'8','green',symbols[2],'white'])
#calc_systems.append([1.129 , 2.623 , 4.530,1.752,r'mrGO-hydroxy-toy model',tcolor,0.0,0.08,1.5,'8','green',symbols[2],colors[0]])
#calc_systems.append([0.919 , 2.447 , 4.416  ,1.752,r'mrGO-hydroxy-edge-model',tcolor,0.0,0.08,1.5,'8','red',symbols[2],colors[1]])
#main result:
#calc_systems.append([0.832 , 2.222 , 4.092 ,1.752,r'mrGO-hydroxy/CeO2(100)',tcolor,0.0,0.08,1.5,'8','red',symbols[2],colors[4]])

#calc_systems.append([1.202 , 2.082 , 4.635  ,1.752,r'mrGO-epoxy-edge-model',tcolor,0.0,0.08,1.5,'8','red',symbols[0],colors[1]])
#calc_systems.append([0.414 , 2.167 , 3.937   ,1.752,r'mrGO-epoxy-patch-model',tcolor,0.0,0.08,1.5,'8','red',symbols[0],colors[3]])

'''
r1='#D21404'
r2='#800000'
r3='#BC544B'
b1='#1338BE'
b2='#1F456E'
b3='#0492C2'
y1='#FFFD37'
y2='#D2B55B'
y3='#E4D00A'
'''
b1='#1338BE'
b2='#7DF9FF'
b3='#4F97A3'
r1='#FC0B00'
r2='#800000'
r3='#FC8900'
g1='#FFFD37'
g2='#F5F5DC'
g3='#CFFF04'

calc_systems.append([1.456, 2.892, 4.432, 1.752, r'CoCo_M1',tcolor,0.0,0.08,1.5,'8',b1,symbols[3],r1])
calc_systems.append([1.270, 2.993, 4.221, 1.752, r'CoCo_M2',tcolor,0.0,0.08,1.5,'8','orchid',symbols[4],r1])
calc_systems.append([1.455, 2.860, 4.552, 1.752, r'CoCu_M1',tcolor,0.0,0.08,1.5,'8','olive',symbols[3],r2])
calc_systems.append([0.829, 2.283, 4.197, 1.752, r'CuCo_M2',tcolor,0.0,0.08,1.5,'8','blue',symbols[2],r2])
calc_systems.append([1.130, 2.460, 4.233, 1.752, r'CoNi_M1',tcolor,0.0,0.08,1.5,'8','pink',symbols[3],r3])
calc_systems.append([1.096, 2.910, 3.962, 1.752, r'NiCo_M2',tcolor,0.0,0.08,1.5,'8','cyan',symbols[2],r3])

calc_systems.append([1.819, 4.105, 4.836, 1.752, r'CuCu_M1',tcolor,0.0,0.08,1.5,'8','black',symbols[3],b1])
calc_systems.append([1.636, 3.981, 4.633, 1.752, r'CuCu_M2',tcolor,0.0,0.08,1.5,'8','black',symbols[2],b1])
calc_systems.append([1.381, 3.642, 4.812, 1.752, r'CuCo_M1',tcolor,0.0,0.08,1.5,'8','blue',symbols[4],b2])
calc_systems.append([1.736, 4.153, 4.695, 1.752, r'CoCu_M2',tcolor,0.0,0.08,1.5,'8','olive',symbols[4],b2])
calc_systems.append([1.817, 4.082, 4.825, 1.752, r'CuNi_M1',tcolor,0.0,0.08,1.5,'8','purple',symbols[4],b3])
calc_systems.append([1.529, 4.235, 4.584, 1.752, r'NiCu_M2',tcolor,0.0,0.08,1.5,'8','green',symbols[2],b3])

calc_systems.append([1.948, 3.941, 4.899, 1.752, r'NiNi_M1',tcolor,0.0,0.08,1.5,'8','orange',symbols[3],g1])
calc_systems.append([1.900, 3.236, 4.742, 1.752, r'NiNi_M2',tcolor,0.0,0.08,1.5,'8','orange',symbols[2],g1])
calc_systems.append([1.944, 3.917, 5.132, 1.752, r'NiCo_M1',tcolor,0.0,0.08,1.5,'8','cyan',symbols[4],g2])
calc_systems.append([1.621, 2.909, 4.441, 1.752, r'CoNi_M2',tcolor,0.0,0.08,1.5,'8','pink',symbols[4],g2])
calc_systems.append([1.764, 3.918, 4.914, 1.752, r'NiCu_M1',tcolor,0.0,0.08,1.5,'8','green',symbols[3],g3])
calc_systems.append([1.971, 3.377, 4.805, 1.752, r'CuNi_M2',tcolor,0.0,0.08,1.5,'8','purple',symbols[2],g3])

#calc_systems.append([0.4605, 0.8505, 3.7208, 1.752, r'Pt(111)',tcolor,0.0,0.08,1.5,'8','orangered',symbols[0],'white'])
'''
import csv
with open('Final_dGs_for_MOOHx_system_ORR.csv', 'w', newline = '') as myfile:
    fieldnames = ['Surface name','dOH','dO','dOOH', 'overpotential', 'onset potential', 'PLS']
    spamwriter = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',', quotechar = '|')
    writer=csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print ('Final dGs for MOOHx_system  system')

    solv_corr_OH=-0.116
    solv_corr_O=-0.083
    solv_corr_OOH=-0.327
    for i in range(len(calc_systems)):
        calc_systems[i][0]+=solv_corr_OH
        calc_systems[i][1]+=solv_corr_O
        calc_systems[i][2]+=solv_corr_OOH
        recalculated_over=overpotential_orr_full(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])
        print  (calc_systems[i][4], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], recalculated_over[0],recalculated_over[1], ' limiting step is ', recalculated_over[2])
        calc_systems[i][3]=recalculated_over[0]
        writer.writerow([calc_systems[i][4], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], recalculated_over[0], recalculated_over[1], recalculated_over[2]])

'''
for i in range(len(calc_systems)):
  if i%2==0 :
     ax.plot(calc_systems[i][0], calc_systems[i][2],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=MarkerStyle('o', fillstyle='left'), markersize=8,label=calc_systems[i][4]+' : %.2f V' % (overpotential_orr(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])))
  
  else :
     ax.plot(calc_systems[i][0], calc_systems[i][2],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=MarkerStyle('o', fillstyle='right'), markersize=8, label=calc_systems[i][4]+' : %.2f V' % (overpotential_orr(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])))

#scaling line:0.8503 x + 3.168
#ax.plot(x,x+3.2,'--',lw=1, dashes=(3,1),c='black')
#ax.plot(x,0.8181*x+3.21,'--',lw=1, dashes=(3,1),c='black')
ax.plot(x,0.87*x+3.22,'--',lw=1, dashes=(3,1),c='black')
ax.text(1.2,2.5 , r'$\Delta$G$_{\sf OOH}$=0.87$\Delta$G$_{\sf OH}$',{'color': 'black', 'fontsize': 10})
ax.text(1.8,2.34 , r'+3.22 eV',{'color': 'black', 'fontsize': 10})
corners=[[1.5,1.2],[x1+(x2-x2)*0.2,y1+(y2-y1)*0.9],[x1+(x2-x2)*0.8,y1+(y2-y1)*0.1],[-2,0]]

#for i in range(len(corners)):
#   ax.text(corners[i][0],corners[i][1], overpotential_label(corners[i][0],corners[i][1]), color='white',fontsize='x-small',horizontalalignment='center',rotation=0,zorder=3)
   
ax.legend(bbox_to_anchor=(-0.15, 1.65), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)

fig.show()
#fig.savefig('ORR_contour_plot_v13_full_scaling_trilayer_v3.pdf', bbox_inches='tight')
# clear the plot
fig.clf()

fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1=-0.7
x2=2.5
ax.axis([x1, x2, x1, ooh_oh_scaling(x2)])


ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$,$\Delta$G$_{\sf O}$ (eV)')

xdata=[]
ydata=[]
y2data=[]

for i in range(len(calc_systems)):
#for i in range(3):
   xdata.append(calc_systems[i][0])
   ydata.append(calc_systems[i][2])
   y2data.append(calc_systems[i][1])


#print #xdata
#print ydata

from pylab import *
fit = polyfit(xdata,ydata,1)
fit_fn = poly1d(fit)
#print fit_fn
aa=fit_fn[1]
bb=fit_fn[0]

fit1 = polyfit(xdata,y2data,1)
fit_fn1 = poly1d(fit1)
#print fit_fn1
#print fit_fn[0], fit_fn[1]
#how bad is scaling
for i in range(len(calc_systems)):
        error=calc_systems[i][2]-(fit_fn[1]*calc_systems[i][0]+fit_fn[0])
        print (error, calc_systems[i])

xx = np.arange(x1,x2, delta)
ax.plot(xx,fit_fn[1]*xx+fit_fn[0],'--',lw=1, dashes=(3,1),c='black',label='OOH scaling')
ax.plot(xx,xx+3.2,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,xx,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,fit_fn1[1]*xx+fit_fn1[0],'--',lw=1, dashes=(3,1),c='red',label='O scaling')


for i in range(len(calc_systems)):
	ax.plot(calc_systems[i][0],calc_systems[i][2],'ro',
           ms=3,marker='o',
           #alpha=0.2,
           color='black')
	ax.plot(calc_systems[i][0],calc_systems[i][1],'ro',
           ms=3,marker='o',
           #alpha=0.2,
           color='red')
	#ax.plot(calc_systems[i][0],calc_systems[i][0],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,label=calc_systems[i][4]+' : %.2f V' % (calc_systems[i][3]),marker=calc_systems[i][11])
	ax.plot(calc_systems[i][0],calc_systems[i][0],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,label=calc_systems[i][4]+' : %.2f V' % (calc_systems[i][3]),marker=calc_systems[i][11])
   	#ax.text(calc_systems[i][0], calc_systems[i][0]+calc_systems[i][7]+0.08,
         #  calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]), color='black',fontsize=6,horizontalalignment='center',
          # rotation=0,zorder=1)

#ax.legend(bbox_to_anchor=(-0.15, 1.55), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)
eq1=r'$\Delta$G$_{\sf OOH}$=%.2f$\Delta$G$_{\sf OH}$+%.2f' %(fit_fn[1],fit_fn[0])
ax.text(-0.5,5.1,eq1,fontsize='small',color='black')
eq2=r'$\Delta$G$_{\sf OOH}$=$\Delta$G$_{\sf OH}$+3.2'
ax.text(-0.5,4.7,eq2,fontsize='small',color='grey')
eq3=r'$\Delta$G$_{\sf O}$=%.2f$\Delta$G$_{\sf OH}$+%.2f' %(fit_fn1[1],fit_fn1[0])
ax.text(-0.5,4.3,eq3,fontsize='small',color='orangered')

#fig.show()
fig.savefig('ORR_scaling_full_scaling_trilayer_v1.pdf', bbox_inches='tight')
# clear the plot
fig.clf()


fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
#fig.patch.set_alpha(0.5)
#fig.patch.set_facecolor('white')
ax.grid(False)
#sns.set_style("white")


#x1=1.23-1
#x2=1.23+1
#y2=
#y1=0

x1=-0.8
x2=2.85
#y2=1
y2=1.23
y1=-0.75

ax.axis([x1, x2, y1, y2])
delta = 0.01
x = np.arange(x1,x2, delta)

ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf O}$ (eV)')

#ax.set_ylabel(r'$\eta_{\sf ORR}$ (V)')
ax.set_ylabel(r'U$_{\sf L, ORR}$ (V)')

#ax.set_ylim(ax.get_ylim()[::-1])


#plot(x,np.maximum(x,1.72-x)-0.5,'--',color='pink',lw=0.67, dashes=(3,1),zorder=2)
#standard scaling
#plot(x,np.minimum(x,-(x+3.2)+4.92),'--',color='grey',lw=0.67, dashes=(3,1),zorder=2)
#calculated scaling
plot(x,np.minimum(x,-(0.87*x + 3.22)+4.92),'--',color='black',lw=0.67, dashes=(3,1),zorder=2)
#plot(x,np.maximum(x,1.72-x)-0.5,'--',color='pink',lw=0.67, dashes=(3,1),zorder=2)
plot(x,np.minimum(x,0.7),'--',color='black',lw=0.67, dashes=(3,1),zorder=2)


#xy=np.array([xp for xp in x if 0.8<xp<0.91])
#ax.fill_between(xy, y2, np.maximum(xy,1.72-xy)-0.5, zorder=1, color='red', alpha=0.3, edgecolor="none")

#for b in x:
#	if(np.maximum(b,3.2-b)-1.23<0.44):
#		print b

#plot(x,np.maximum(x,bb-x*(aa)-0.65)-1.23,'--',color='grey',lw=0.67, dashes=(3,1))
#slope not one
#plot(x,np.maximum(x,3.18-0.82*x-0.35)-1.23,'--',color='pink',lw=0.67,dashes=(3,1))
#plot(x,np.maximum(x,2.46-x)-1.23,'-',color='black',lw=0.67)

#import matplotlib.patches as patches
#ax.add_patch(
#    patches.Rectangle(
#        (calc_systems[1][1]-calc_systems[1][0],calc_systems[1][3]-0.04),   # (x,y)
#        0.25,          # width
#        0.05,          # height
#        fill=True,
#        edgecolor="none", 
#        facecolor='red',
#    )
#)


   #ax.plot(calc_systems[i][0],1.23-calc_systems[i][3],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=calc_systems[i][11],label=calc_systems[i][4]+' : %.2f V' % (1.23-calc_systems[i][3]))

for i in range(len(calc_systems)):
   over_cal=overpotential_orr(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])

   if i%2==0 :
      ax.plot(calc_systems[i][0], 1.23-over_cal,calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=MarkerStyle('D', fillstyle='full'), markersize=6,label=calc_systems[i][4]+' : %.2f V' % (1.23-over_cal))

   else :
      ax.plot(calc_systems[i][0], 1.23-over_cal, calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=MarkerStyle('s', fillstyle='full'), markersize=6,label=calc_systems[i][4]+' : %.2f V' % (1.23-over_cal))


   #ax.plot(calc_systems[i][0],1.23-over_cal,calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker='o',label=calc_systems[i][4]+' : %.2f V' % (1.23-over_cal))
   #if(i!=1):
   #	ax.text(calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
#		color='black',fontsize=6,horizontalalignment='left',rotation=0,zorder=4)
 #  else:
  # 	ax.text(calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
#		color='black',fontsize=6,horizontalalignment='right',rotation=0,zorder=4)
ax.legend(bbox_to_anchor=(-0.15, 1.4255), loc=2, borderaxespad=0, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)
fig.savefig('ORR_1D_plot_fix_trilayer_v3.pdf', bbox_inches='tight')
fig.show('ORR_1D_plot_full_scaling_trilayer_v3.pdf')

for i in range(len(calc_systems)):
	print ('%s	a	OH	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][0]))
	print ('%s	a	O	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][1]))
	print ('%s	a	OOH	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][2]))
for i in range(len(calc_systems)):	
	print ("'%s'," %(calc_systems[i][4]))

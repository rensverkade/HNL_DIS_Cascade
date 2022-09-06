#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py as h5

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT10nlo/0")
p = lhapdf.mkPDF("CT18NLO/0")

s2_W =0.22290 #sin_W^2
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g

hf = h5.File('Struct_func.hdf5', 'r')
xs_x = np.array(hf.get('xs_x'))
qs_x =  np.array(hf.get('qs_x'))
xs_Q = np.array(hf.get('xs_Q'))
qs_Q =  np.array(hf.get('qs_Q'))
F2_x_y =  np.array(hf.get('F2_x_y'))
F2_x_yZ = np.array( hf.get('F2_x_yZ'))
F2_x_Z =  np.array(hf.get('F2_x_Z'))
F3_x_y = np.array( hf.get('F3_x_y'))
F3_x_yZ =  np.array(hf.get('F3_x_yZ'))
F3_x_Z = np.array( hf.get('F3_x_Z'))

F2_Q_y =  np.array(hf.get('F2_Q_y'))
F2_Q_yZ = np.array( hf.get('F2_Q_yZ'))
F2_Q_Z =  np.array(hf.get('F2_Q_Z'))
F3_Q_y = np.array( hf.get('F3_Q_y'))
F3_Q_yZ =  np.array(hf.get('F3_Q_yZ'))
F3_Q_Z = np.array( hf.get('F3_Q_Z'))



hf.close()
l=len(xs_Q)
print(l)
k = len(xs_x)
close = True
#close = False
colours = np.array(['r','b','g','orange','purple','black','pink','cyan'])
#e_q = np.array([-1/3,2/3,-1/3,2/3,-1/3])
e_q = np.array([1/3,2/3,1/3,2/3])


#----------------
############
# Time for #
# F1 F2 FU #
############


G_F = 1.188*10**(-5) #GeV^-2
M_W =  80.38 #GeV
mp = 0.94 #GeV
g_A = np.array([-0.5, 0.5,-0.5,0.5])
g_V = g_A -2*e_q*s2_W
fig, axs = plt.subplots(2, 2, sharex=True)
for i in range(l): #l 0-8        
    axs[0,0].plot(xs_x,F2_x_y[i,:],color=colours[i], label = r'$Q^2=$'+str(qs_x[i]))
    axs[1,0].plot(xs_x,F2_x_yZ[i,:],color=colours[i], ls='dashed')
    axs[0,1].plot(xs_x,F2_x_Z[i,:],color=colours[i], ls= 'dotted')
    if  i==5:     
        axs[1,1].plot(xs_x,F2_x_y[i,:],color=colours[i], label = r'$F^{\gamma}_2$')
        axs[1,1].plot(xs_x,F2_x_yZ[i,:],color=colours[i], label = r'$F^{\gamma Z}_2$', ls='dashed')
        axs[1,1].plot(xs_x,F2_x_Z[i,:],color=colours[i], label = r'$F_2^Z$', ls= 'dotted')

axs[1,0].set_xlabel(r'$x$')
axs[1,1].set_xlabel(r'$x$')
axs[0,0].set_ylabel(r'$F2 (x, Q^2)$')
axs[1,0].set_ylabel(r'$F2 (x, Q^2)$')

axs[0,0].set_ylim(-0.1,3)
axs[1,0].set_ylim(-0.1,5)
axs[0,1].set_ylim(ymin=-1)
#axs[1,1].set_ylim(r'$F2 (x, Q^2)$')
plt.xscale('log')
#plt.yscale('log')
#plt.ylim(0,1.5)
fig.tight_layout()
fig.legend(loc='center right')
plt.subplots_adjust(hspace=0.1,right=0.75)
plt.savefig('Struct_figs_SSM/F2_p_x.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()
#plt.show()

fig, axs = plt.subplots(2, 2, sharex=True)
for i in range(l): #l 0-8        
    axs[0,0].plot(qs_Q,F2_Q_y[i,:],color=colours[i], label = r'$x=$'+str(xs_Q[i]))
    axs[1,0].plot(qs_Q,F2_Q_yZ[i,:],color=colours[i], ls='dashed')
    axs[0,1].plot(qs_Q,F2_Q_Z[i,:],color=colours[i], ls= 'dotted')
    if  i==5:     
        axs[1,1].plot(qs_Q,F2_Q_y[i,:],color=colours[i], label = r'$F^{\gamma}_2$')
        axs[1,1].plot(qs_Q,F2_Q_yZ[i,:],color=colours[i], label = r'$F^{\gamma Z}_2$', ls='dashed')
        axs[1,1].plot(qs_Q,F2_Q_Z[i,:],color=colours[i], label = r'$F_2^Z$', ls= 'dotted')

axs[1,0].set_xlabel(r'$Q^2$')
axs[1,1].set_xlabel(r'$Q^2$')
axs[0,0].set_ylabel(r'$F2 (x, Q^2)$')
axs[1,0].set_ylabel(r'$F2 (x, Q^2)$')
plt.xscale('log')
axs[0,0].set_yscale('log')
axs[0,1].set_yscale('symlog')
axs[1,0].set_yscale('log')
axs[1,1].set_yscale('log')
#plt.ylim(0,1.5))
fig.tight_layout()
fig.legend(loc='center right')
plt.subplots_adjust(hspace=0.1,right=0.75)
plt.savefig('Struct_figs_SSM/F2_p_Q2.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

    
quit()




#=========

##############
# Time to add#
#       Z   #
##############










## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

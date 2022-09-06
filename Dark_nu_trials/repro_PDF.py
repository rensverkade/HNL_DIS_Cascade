#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
p = lhapdf.mkPDF("CT10nlo", 0)
p = lhapdf.mkPDF("CT10nlo/0")


## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g


# TODO: demonstrate looping over PDF set members
#pset = lhapdf.getPDFSet("CT10nlo")
#print('descript')
#print(pset.description)
#pcentral = pset.mkPDF(0)
#pdfs1 = pset.mkPDFs()
#pdfs2 = lhapdf.mkPDFs("CT10nlo") # a direct way to get all the set's PDFs

## Filling a numpy 2D array with xf(x,Q) samples

xs_x = [x for x in np.logspace(-4, 0, 80)]
qs_Q = [q for q in np.logspace(0, 4, 80)]
xs_Q= np.array([10**-8,10**-6, 0.0001,0.01,0.1,0.2,0.5,0.8])
qs_x = np.array([10,50,100,200,500,1000,2000,5000])
gluon_xfs_x = np.empty([len(qs_x), len(xs_x)])
#up_xfs_x = np.empty([len(xs), len(qs)])
#down_xfs_x = np.empty([len(xs), len(qs)])
gluon_xfs_Q = np.empty([len(xs_Q), len(qs_Q)]) #(8,80)
#up_xfs_Q = np.empty([len(xs), len(qs)])
#down_xfs_Q = np.empty([len(xs), len(qs)])
l=len(xs_Q)
print(l)
k = len(xs_x)
close = True
colours = np.array(['r','b','g','orange','purple','pink','black','cyan'])


for i in range(l): #0-8
    for j in range(k): #0-80   p.xfxQ2 takes Q^2 as arg
        gluon_xfs_Q[i,j] = p.xfxQ(21, xs_Q[i], qs_Q[j])#taking Q^2 as argument
        gluon_xfs_x[i,j] = p.xfxQ(21, xs_x[j], qs_x[i])#taking Q^2 as argument
        #up_xfs[ix,iq] = p.xfxQ2(1, x, q)#taking Q^2 as argument
        #down_xfs[ix,iq] = p.xfxQ2(2, x, q)#taking Q^2 as argument
    plt.plot(xs_x,gluon_xfs_x[i,:],color=colours[i], label = r'flav = g, $Q^2=$'+str(qs_x[i]))

plt.xlabel(r'$x$')
plt.ylabel(r'$xf(x,Q)$')
plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.savefig('PDF_figs_repro/PDF_g_x.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

for i in range(l):
    plt.plot(qs_Q,gluon_xfs_Q[i,:],color=colours[i], label = r'flav = g, $x=$'+str(xs_Q[i]))

plt.xlabel(r'$Q$')
plt.ylabel(r'$xf(x,Q)$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('PDF_figs_repro/PDF_g_Q.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()


    


#plt.close()


## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

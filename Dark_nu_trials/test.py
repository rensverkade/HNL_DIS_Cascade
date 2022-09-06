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

## Gluon PDF querying at x=0.001, Q2=10000
print('gluon,  x=0.001, Q2=10000')
print(p.xfxQ2(21, 1e-3, 1e4))

## Basic all-flavour PDF querying at x=0.01, Q=M_Z
print('flavours')
for pid in p.flavors():
    print(str(pid))
    print(p.xfxQ(pid, 0.01, 91.2))
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g


# TODO: demonstrate looping over PDF set members
pset = lhapdf.getPDFSet("CT10nlo")
print('descript')
print(pset.description)
pcentral = pset.mkPDF(0)
pdfs1 = pset.mkPDFs()
pdfs2 = lhapdf.mkPDFs("CT10nlo") # a direct way to get all the set's PDFs

## Filling a numpy 2D array with xf(x,Q) samples

xs = [x for x in np.logspace(-7, 0, 40)]
qs = [q for q in np.logspace(1, 16, 50)]
gluon_xfs = np.empty([len(xs), len(qs)])
up_xfs = np.empty([len(xs), len(qs)])
down_xfs = np.empty([len(xs), len(qs)])
for ix, x in enumerate(xs):
    for iq, q in enumerate(qs):
        gluon_xfs[ix,iq] = p.xfxQ2(21, x, q)#taking Q^2 as argument
        up_xfs[ix,iq] = p.xfxQ2(1, x, q)#taking Q^2 as argument
        down_xfs[ix,iq] = p.xfxQ2(2, x, q)#taking Q^2 as argument
print(gluon_xfs)

close = True

plt.imshow(gluon_xfs, origin='lower', cmap=cm.hot)#,interpolation='bicubic')
plt.colorbar(label=r'$xf(x,Q^2$)')
plt.xlabel('x')
plt.ylabel(r'$Q^2$')
#plt.xscale('log')
#plt.yscale('log')
plt.title('PDF vals gluon?')
plt.savefig('PDF_figs_test/test_PDF_gluon_2d.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

plt.imshow(up_xfs, origin='lower', cmap=cm.hot)#,interpolation='bicubic')
plt.colorbar(label=r'$xf(x,Q^2$)')
plt.xlabel('x')
plt.ylabel(r'$Q^2$')
#plt.xscale('log')
#plt.yscale('log')
plt.title('PDF vals up?')
plt.savefig('PDF_figs_test/test_PDF_up_2d.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

plt.imshow(down_xfs, origin='lower', cmap=cm.hot)#,interpolation='bicubic')
plt.colorbar(label=r'$xf(x,Q^2$)')
plt.xlabel('x')
plt.ylabel(r'$Q^2$')
#plt.xscale('log')
#plt.yscale('log')
plt.title('PDF vals down?')
plt.savefig('PDF_figs_test/test_PDF_down_2d.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

plt.plot(gluon_xfs[30,:])
plt.xlabel(r'$Q^2$')
plt.ylabel('PDF val')
plt.title('PDF vals gluon for x='+str(xs[3]))
plt.xscale('log')
plt.yscale('log')
plt.savefig('PDF_figs_test/test_PDF_gluon_1d_Q.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

plt.plot(up_xfs[30,:])
plt.xlabel(r'$Q^2$')
plt.ylabel('PDF val')
plt.title('PDF vals u for x='+str(xs[3]))
plt.xscale('log')
#plt.yscale('log')
plt.savefig('PDF_figs_test/test_PDF_up_1d_Q.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()


plt.plot(down_xfs[30,:])
plt.xlabel(r'$Q^2$')
plt.ylabel('PDF val')
plt.title('PDF vals d for x='+str(xs[3]))
plt.xscale('log')
#plt.yscale('log')
plt.savefig('PDF_figs_test/test_PDF_down_1d_Q.pdf', dpi=800,bbox_inches = 'tight')
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

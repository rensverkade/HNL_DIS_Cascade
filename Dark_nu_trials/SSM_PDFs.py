#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT10nlo/0")
p = lhapdf.mkPDF("CT18NLO/0")
pw = lhapdf.mkPDF("CT10wnlo/0")


## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g



xs_x = [x for x in np.logspace(-4, 0, 80)]
qs_Q = [q for q in np.logspace(0, 4, 80)]
xs_Q= np.array([10**-8,10**-6,5*10**-6,10**-5, 0.0001,0.01,0.1,0.5])
qs_x = np.array([10,50,100,200,500,1000,2000,5000])
qs_x2 = np.array([0.05,0.1,1,5,10,50,100,200,500])
xs_Q2= np.array([10**-10,10**-9,10**-8,10**-6, 0.0001,0.01,0.1])


gluon_xfs_x = np.zeros([len(qs_x), len(xs_x)])
up_xfs_x = np.zeros([len(qs_x), len(xs_x)])
down_xfs_x = np.zeros([len(qs_x), len(xs_x)])
gluon_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
up_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)])
down_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)])


up_b_xfs_x = np.zeros([len(qs_x), len(xs_x)])
down_b_xfs_x = np.zeros([len(qs_x), len(xs_x)])
up_b_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)])
down_b_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)])


l=len(xs_Q)
print(l)
k = len(xs_x)
close = True
colours = np.array(['r','b','g','orange','purple','black','pink','cyan'])
#e_q = np.array([-1/3,2/3,-1/3,2/3,-1/3])
e_q = np.array([1/3,2/3,1/3,2/3,1/3])
#e_q = np.array([2/3,1/3,1/3,2/3,1/3])

fig, axs = plt.subplots(2, 2, sharex=True)
for i in range(l): #0-8
    for j in range(k): #0-80   p.xfxQ2 takes Q^2 as arg
        #gluon_xfs_Q[i,j] = p.xfxQ2(21, xs_Q2[i], qs_Q[j])#taking Q^2 as argument
        #gluon_xfs_x[i,j] = p.xfxQ2(21, xs_x[j], qs_x2[i])#taking Q^2 as argument
        up_xfs_x[i,j] = p.xfxQ2(2, xs_x[j], qs_x[i])#taking Q^2 as argument
        down_xfs_x[i,j] = p.xfxQ2(1, xs_x[j], qs_x[i])#taking Q^2 as argument
        up_xfs_Q[i,j] = p.xfxQ2(2, xs_Q[i], qs_Q[j])#taking Q^2 as argument
        down_xfs_Q[i,j] = p.xfxQ2(1, xs_Q[i], qs_Q[j])#taking Q^2 as argument
        up_b_xfs_x[i,j] = p.xfxQ2(-2, xs_x[j], qs_x[i])#taking Q^2 as argument
        down_b_xfs_x[i,j] = p.xfxQ2(-1, xs_x[j], qs_x[i])#taking Q^2 as argument
        up_b_xfs_Q[i,j] = p.xfxQ2(-2, xs_Q[i], qs_Q[j])#taking Q^2 as argument
        down_b_xfs_Q[i,j] = p.xfxQ2(-1, xs_Q[i], qs_Q[j])#taking Q^2 as argument
for i in range(1,l,2): #0-8
    axs[0,0].plot(xs_x,up_xfs_x[i,:],color=colours[i], label = r'$Q^2=$'+str(qs_x[i]))
    axs[0,1].plot(xs_x,down_xfs_x[i,:],color=colours[i],ls='dashed')
    axs[0,0].plot(xs_x,up_b_xfs_x[i,:],color=colours[i],ls='dotted')
    axs[0,1].plot(xs_x,down_b_xfs_x[i,:],color=colours[i],ls='dashdot')
    if i==5:
        axs[1,0].plot(xs_x,up_xfs_x[i,:],color=colours[i], label = 'flav = u')
        axs[1,0].plot(xs_x,down_xfs_x[i,:],color=colours[i],ls='dashed',label='flav = d')
        axs[1,1].plot(xs_x,up_b_xfs_x[i,:],color=colours[i],ls='dotted', label = r'flav = $\overline{u}$')
        axs[1,1].plot(xs_x,down_b_xfs_x[i,:],color=colours[i],ls='dashdot',label=r'flav = $\overline{d}$')
        
axs[1,1].set_xlabel(r'$x$')
axs[1,0].set_xlabel(r'$x$')
axs[0,0].set_ylabel(r'$xf(x,Q^2)$')
axs[1,0].set_ylabel(r'$xf(x,Q^2)$')
plt.xscale('log')
#plt.yscale('log')
fig.tight_layout()
fig.legend(loc='center right')
plt.subplots_adjust(hspace=0.1,right=0.75)
plt.savefig('PDF_figs_SSM/PDFs_x.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()


fig, axs = plt.subplots(2, 2, sharex=True)
for i in range(1,l,2):
    axs[0,0].plot(qs_Q,up_xfs_Q[i,:],color=colours[i], label = r' $x=$'+str(xs_Q[i]))
    axs[0,1].plot(qs_Q,down_xfs_Q[i,:],color=colours[i],ls='dashed')
    axs[0,0].plot(qs_Q,up_b_xfs_Q[i,:],color=colours[i],ls='dotted')
    axs[0,1].plot(qs_Q,down_b_xfs_Q[i,:],color=colours[i],ls='dashdot')

    if i==5:
        axs[1,0].plot(qs_Q,up_xfs_Q[i,:],color=colours[i], label = 'flav = u')
        axs[1,0].plot(qs_Q,down_xfs_Q[i,:],color=colours[i],ls='dashed',label='flav = d')
        axs[1,1].plot(qs_Q,up_b_xfs_Q[i,:],color=colours[i],ls='dotted', label = r'flav = $\overline{u}$')
        axs[1,1].plot(qs_Q,down_b_xfs_Q[i,:],color=colours[i],ls='dashdot',label=r'flav = $\overline{d}$')
    

axs[1,0].set_xlabel(r'$Q^2$')
axs[1,1].set_xlabel(r'$Q^2$')
axs[0,0].set_ylabel(r'$xf(x,Q^2)$')
axs[1,0].set_ylabel(r'$xf(x,Q^2)$')
plt.xscale('log')
axs[0,0].set_yscale('symlog')
axs[0,1].set_yscale('symlog')
#axs[1,0].set_yscale('symlog')
#axs[1,1].set_yscale('symlog')
fig.tight_layout()
fig.legend(loc='center right')
plt.subplots_adjust(hspace=0.1,right=0.75)
plt.savefig('PDF_figs_SSM/PDFs_Q.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

quit()

#----------------
############
# Time for #
# F1 F2 FU #
############


G_F = 1.188*10**(-5) #GeV^-2
M_W =  80.38 #GeV
mp = 0.94 #GeV 


def struct_PDF(x,Q2):
    flavs = p.flavors()[:-1]
    
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(p.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           anti-quarks
    xf_g = np.array(p.xfxQ2(21, x, Q2))
    #print('xf',xf)
    #F2 = np.zeros(f_len)
    #F2 = np.sum(xf[f_half:]+xf[:f_half]) #sea qu
    
    F2 = np.sum(e_q**2 *(xf[f_half:]+xf[:f_half])) #sea quarks+valence, ?c dominated??
    
    xF3 = np.sum((xf[f_half:]-xf[:f_half])) #valence quarks(sea cancels)#times 2 for neutrino
    FL=0 #WHERE am is supposeed to find disnshit
    x2F1 = F2*(1+4*mp**2 * x**2 /Q2)/(1+FL)
    return (x2F1,F2,xF3)

'''
def struct_PDFw(x,Q2):
    flavs = pw.flavors()[:-1]
    
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(pw.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           anti-quarks
    xf_g = np.array(pw.xfxQ2(21, x, Q2))
    #print('xf',xf)
    F2 = np.zeros(f_len)
    F2 = 2*np.sum(xf[f_half:]+xf[:f_half]) #sea quarks+valence, ?c dominated??
    xF3 = 2*np.sum(xf[f_half:]-xf[:f_half]) #valence quarks(sea cancels)
    FL=0 #WHERE am is supposeed to find disnshit
    x2F1 = F2*(1+4*mp**2 * x**2 /Q2)/(1+FL)
    return (x2F1,F2,xF3)
'''
#qs_x = np.array([3.5,90])

F1_x = np.zeros([len(qs_x), len(xs_x)])
F2_x = np.zeros([len(qs_x), len(xs_x)])
F3_x = np.zeros([len(qs_x), len(xs_x)])
F1_Q = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F2_Q = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F3_Q = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)


for i in range(l): #l 0-8
    for j in range(k): #0-80   p.xfxQ2 takes Q^2 as arg
        F1_Q[i,j], F2_Q[i,j], F3_Q[i,j] = struct_PDF(xs_Q[i], qs_Q[j])#taking Q^2 as argument
        F1_x[i,j], F2_x[i,j] ,F3_x[i,j]  = struct_PDF(xs_x[j], qs_x[i])#taking Q^2 as argument
        
    plt.plot(xs_x,F2_x[i,:],color=colours[i], label = r'F2, $Q^2=$'+str(qs_x[i]))

plt.xlabel(r'$x$')
plt.ylabel(r'$F2 (x, Q^2)$')
plt.xscale('log')
#plt.yscale('log')
#plt.ylim(0,1.5)
plt.legend()
plt.savefig('PDF_figs_SSM/F2_x.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()

for i in range(l):
    plt.plot(qs_Q,F2_Q[i,:],color=colours[i], label = r'F2, $x=$'+str(xs_Q[i]))

plt.xlabel(r'$Q$')
plt.ylabel(r'$F2(x,Q^2)$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('PDF_figs_SSM/F2_Q2.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()


#=========

##############
# Time to add#
#       Z'   #
##############










## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

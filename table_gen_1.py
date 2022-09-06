import numpy as np
import matplotlib.pyplot as plt
import math
import h5py as h5
import sys
from numpy import linalg as la
hf = h5.File('HNL_mass_pascoli.hdf5','r')

m_D = np.array(hf.get('m_D'))
L_L = np.array(hf.get('L_L'))
L_R = np.array(hf.get('L_R'))
mu= np.array(hf.get('mu'))
M_X = np.array(hf.get('M_X'))
M_diag = np.array(hf.get('M_diag'))
M_mix = np.array(hf.get('mix'))

hf.close()
print(L_L)
print(np.shape(M_diag))

def npy_to_tex(name, A, var_name=None):
    """Save to LaTeX tabular format"""
    if not var_name:
        var_name = name
    A = A.tolist()
    tex_string = [r"\begin{table} \n\centering \n \begin{table}"]
    for col in A:
        tex_string += "|c"
    tex_string += "|"
    
    text_list = []
    for row in A:
        tex_elements = []
        if type(row)is not list:
            text_list.append(r"{:.2e}".format(row))
        else: 
            for element in row:
                tex_elements.append("{:.2e}".format(element))
            text_list.append(" & ".join(tex_elements))
    tex_string.append("\\\\ \n".join(text_list))
    tex_string += "\end{tabular}"
    tex_string += "\caption{}"
    tex_string += "\end{table}"
    text_string = "".join(tex_string)
    file  = open("values/{}.tex".format(name), "wt")
    file.write(text_string)
    file.close()
    
npy_to_tex('\Lambda_L',L_L)
npy_to_tex('\Lambda_R',L_R)
npy_to_tex('mu',mu)
npy_to_tex('M',M_diag)
npy_to_tex('M_X',M_X)

for i in range(4):
    npy_to_tex('M_mix_'+str(i),M_mix[i,:,:]*M_mix[i,:,:])

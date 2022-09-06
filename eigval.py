import numpy as np
import matplotlib.pyplot as plt
import math
import h5py as h5
import sys
from numpy import linalg as la
e_m = 0.001


def Brent(xa,xb,f,arg1,arg2):
    """Returns a root of input function f(x) within target error using Brent's Method.
       range: [xa,xb]
       
       Adapted from Wikipedia pseudocode:
       https://en.wikipedia.org/wiki/Brent%27s_method#Algorithm
    """
    print("Finding Root (Brent) ...")
    i = 0
    tgt_error=0.0001
    mflag = True
    xc = xa
    while abs(xa-xb) > tgt_error:
        fa = f(xa,arg1,arg2)
        fb = f(xb,arg1,arg2)
        fc = f(xc,arg1,arg2)
        if fa != fc and fb != fc:
            s = (xa*fb*fc)/((fa-fb)*(fa-fc)) + (xb*fa*fc)/((fb-fa)*(fb-fc)) \
              + (xc*fa*fb)/((fc-fa)*(fc-fb))
        else:
            s = xb - fb*(xb-xa)/(fb-fa)
        if ((s < (3*xa+xb)*0.25 or s > xb) or \
            (mflag == True and abs(s-xb) >= abs(xb-xc)*0.5) or \
            (mflag == False and abs(s-xb) >= abs(xc-xd)*0.5) or \
            (mflag == True and abs(xb-xc) < abs(tgt_error)) or \
            (mflag == False and abs(xc-xd) < abs(tgt_error))):
            s = 0.5*(xa+xb)
            mflag = True
        else:
            mflag = False
        fs = f(s,arg1,arg2)
        xd = xc
        xc = xb
        if fa*fs < 0:
            xb = s
        else:
            xa = s
        if abs(fa) < abs(fb):
            xa,xb = xb,xa
        i += 1
        sys.stdout.write("\r  Steps = "+str(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return xb

def bisection(x_min,x_max,func,arg1,arg2,ref): #bisection root finding method
    a =x_min
    b =x_max
    s1 = func(a,arg1,arg2)*func(b,arg1,arg2)
    step = 0
    roots = np.empty(0)
    i=0
    if s1<0:
        print('there is a root')#did our initial bracket contain 1 root?
    else:
        print('bad start') #if not we chose wrong
        return bisection(x_min,x_max*0.75,func,arg1,arg2,ref)
    while math.sqrt((a-b)*(a-b)) > ref:
        if i == 0:
            c  = (a+b) /2#generate new point halfway between a and b
        else:
            c/=20#adjust
        if func(a,arg1,arg2)*func(c,arg1,arg2)<0: # is the root within [a,c]
            print('b')
            b=c #b is replaced
            i=0
        elif (func(b,arg1,arg2)*func(c,arg1,arg2)<0):# or is the root witgin [c,b]?
            print('a')
            a=c #a is replaced
            i=0
        else:
            print(i)
            i=i+1
            print(i)
            if i > 20:
                print('bisection failed')
                return False,step
        step +=1 #keep track of the nr of steps
    return((a+b)/2),step

def det(A):
    # Section 1: Establish n parameter and copy A
    n = len(A)
    AM =np.copy(A)
 
    # Section 2: Row ops on A to get in upper triangle form
    for fd in range(n): # A) fd stands for focus diagonal
        for i in range(fd+1,n): # B) only use rows below fd row
            #if AM[fd][fd] == 0: # C) if diagonal is zero ...
            #    AM[fd][fd] == 1.0e-18 # change to ~zero
            # D) cr stands for "current row"
            crScaler = AM[i][fd] / AM[fd][fd] 
            # E) cr - crScaler * fdRow, one element at a time
            for j in range(n): 
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
     
    # Section 3: Once AM is in upper triangle form ...
    product = 1.0
    for i in range(n):
        # ... product of diagonals is determinant
        product *= AM[i][i] 
 
    return product

def det_r(A):
    # Section 1: store indices in list for row referencing
    indices = list(range(len(A)))
     
    # Section 2: when at 2x2 submatrices recursive calls end
    if len(A) == 2 and len(A[0]) == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val
 
    # Section 3: define submatrix for focus column and 
    #      call this function
    for fc in indices: # A) for each focus column, ...
        # find the submatrix ...
        As = np.copy(A) # B) make a copy, and ...
        As = As[1:] # ... C) remove the first row
        height = len(As) # D) 
 
        for i in range(height): 
            # E) for each remaining row of submatrix ...
            #     remove the focus column elements
            As[i] = As[i][0:fc] + As[i][fc+1:] 
 
        sign = (-1) ** (fc % 2) # F) 
        # G) pass submatrix recursively
        sub_det = determinant_recursive(As)
        # H) total all returns from recursion
        total += sign * A[0][fc] * sub_det 
 
    return total
  

def char_eq(l,A,r):
    n = np.shape(A)[0]
    B = A-np.diag(np.ones(n)*l)
    if r is None:
        d = la.det(B)
    else:
        print(r)
        #print('prod',np.prod(l-r))
        d = la.det(B)/(np.prod(l-r))
        #print('d',d)
    return d




y_nu = np.ones((3,1))
#print(y_nu[1][0])
print(np.shape(y_nu), np.shape(y_nu.T))
print(y_nu)
N_hnl=2 # nr of HNLS

m_n =1
v_H = 2.46 #10^2 GeV
v_phi = 2.40 

y_nu = np.ones((3,1)).T #times yukawa coupling
y_N =np.ones((N_hnl,1)).T #times yukawa coupling
m_D = np.zeros((N_hnl,3))
m_D[0,:] = v_H/np.sqrt(2)  *y_nu #(N,3) #one Md and filled with zeros
mu =np.ones((N_hnl,N_hnl))*m_n
L = v_phi/np.sqrt(2)*y_N #(1,N)

A = np.zeros((3,3))
#print(m_D,m_D.T,L,L.T, mu)

print(np.shape(m_D.T))
M = np.block([ [np.zeros((3,3)),m_D.T,np.zeros((3,1)) ],
               [m_D,mu,L.T ],
               [np.zeros((1,3)),L,0]])
print(np.round(M))
print(la.eig(M))
quit()

B = np.array([[7,3,0],[3,3,0],[6,3,5]])
#B = np.array([[7,3],[3,3]])

print(B)
ls = np.arange(-10,30,0.5)

l0,steps = bisection(-10,30,char_eq,B,None,0.01) 

print(l0)
chars = np.zeros(len(ls))

for i in range(len(ls)):
    chars[i] = char_eq(ls[i],B,None)
#plt.plot(ls,chars)
#plt.scatter(l0,char_eq(l0,B,None))
#plt.show()

#plt.plot(ls,chars)
#plt.show()

ms = np.zeros(3)
r=np.array([])
'''
for i in range(3):
    ms[i], steps = bisection(-10,30,char_eq,B,r,0.01) 
    r = np.append(r,ms[i])
    plt.scatter(ms[i],char_eq(ms[i],B,None))

plt.plot(ls,chars)
#plt.scatter(ms[0],char_eq(ms[0],B,None))
#plt.scatter(ms[1],char_eq(ms[1],B,None))
plt.show()
'''          

#quit()
N = len(M[0,:])
print('N',N)
ms = np.zeros(N)
r=np.array([])
#N=3

for i in range(N):
    ms[i] = Brent(1,70,char_eq,M,r) 
    plt.scatter(ms[i],char_eq(ms[i],M,None))

    ls = np.arange(-10,100,5)
    chars = np.zeros(len(ls))
    for j in range(len(ls)):
        chars[j] = char_eq(ls[j],M,r)
    
    plt.plot(ls,chars, label=str(i))
    
    r = np.append(r,ms[i])
    
print(ms)
#plt.scatter(ms[0],char_eq(ms[0],B,None))
#plt.scatter(ms[1],char_eq(ms[1],B,None))
plt.legend()
ymin = -0.2*10**15
#plt.ylim(ymin,10**16)

#plt.xlim(-100,500)
#plt.savefig('plts/eig.pdf',dpi=800)
plt.show()

for i in range(N):
    plt.scatter(ms[i],char_eq(ms[i],M,r[0:i]))

    ls = np.arange(-10,60,5)
    chars = np.zeros(len(ls))
    for j in range(len(ls)):
        chars[j] = char_eq(ls[j],M,r[0:i])
    
    plt.plot(ls,chars, label=str(i))
    plt.show()
    r = np.append(r,ms[i])


quit()


ls = np.arange(-100,800,5)
chars = np.zeros(len(ls))
for j in range(len(ls)):
    chars[j] = char_eq(ls[j],M,r[0:3])

plt.scatter(ms[3],char_eq(ms[3],M,r[0:3]))
plt.plot(ls,chars)
plt.show()

print(ms,r[0:3])





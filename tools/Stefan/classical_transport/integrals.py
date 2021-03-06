from __future__ import division, print_function

from scipy.special import erf
from scipy import exp,sqrt
from scipy.constants import pi
from scipy.integrate import quad

def Chandrasehkar(x):
    return (erf(x) - 2*x*exp(-x**2)/sqrt(pi))/(2*x**2)

def fC1(x,y):
    return y * x**3 * Chandrasehkar(x) * exp(-y*x**2)

def fC2(x,y):
    return (y**2 * x**5 - y* x**3) * Chandrasehkar(x) * exp(-y*x**2)

def FC1(y):
    a,err = quad(fC1,0,float("inf"),args=(y,))
    return a

def FC2(y):
    a,err = quad(fC2,0,float("inf"),args=(y,))
    return a

def f(x,y):
    return erf(x) * x * exp(-y*x**2)

def f2(z,y):
    return erf(z/sqrt(y)) * z * exp(-z**2)

def h(x,y):
    return erf(x) * x**3 * exp(-y*x**2)

def h2(z,y):
    return erf(z/sqrt(y)) * z**3 * exp(-z**2)

def k(x,y):
    return (y*x**3 - 3*x/2) * erf(x) * exp(-y*x**2)

def k2(z,y):
    return (z**3 - 3*z/2) * erf(z/sqrt(y)) * exp(-z**2)

def F(y):
    ## Alternative solution: use the fact that erf(x) is ~1 for large x
    ## numerical for [0,b]
    #a,err = quad(f,0,b,args=(y,))
    ## add analytical for [b,inf] for erf(x) ~ 1
    #a = a + exp(-b**2 * y)/(2*y)
    if y<1:
        a,err = quad(f,0,float("inf"),args=(y,))
        a = y*0.5*a
    else:
        a,err = quad(f2,0,float("inf"),args=(y,))
        a = 0.5*a
    return a - y/(4*(1+y)**(3/2))

def H(y):
    a,err = quad(h,0,float("inf"),args=(y,))
    ## Alternative solution: use the fact that erf(x) is ~1 for large x
    ## numerical for [0,b]
    #a,err = quad(f,0,b,args=(y,))
    ## add analytical for [b,inf] for erf(x) ~ 1
    #a = a + exp(-b**2 * y) * (b**2 * y + 1)/(2*y**2)
    if y<1:
        a,err = quad(h,0,float("inf"),args=(y,))
        a = y**2*0.5*a
    else:
        a,err = quad(h2,0,float("inf"),args=(y,))
        a = 0.5*a
    return a - (y**2)*3/(8*(1+y)**(5/2))

def K(y):
    a,err = quad(k,0,float("inf"),args=(y,))
    ## Alternative solution: use the fact that erf(x) is ~1 for large x
    ## numerical for [0,b]
    #a,err = quad(f,0,b,args=(y,))
    ## add analytical for [b,inf] for erf(x) ~ 1
    #a = a + exp(-b**2 * y) * (b**2 * y + 1)/(2*y**2)
    if y<1:
        a,err = quad(k,0,float("inf"),args=(y,))
        a = y*0.5*a
    else:
        a,err = quad(k2,0,float("inf"),args=(y,))
        a = 0.5*a
    return a + y*3/(8*(1+y)**(5/2))

def Fa(y):
    #analytical form of F(y), derived by H.M. Smith 2018-04-20
    return 1/(4*(1+y)**(3/2))

def Ha(y):
    #analytical form of H(y), derived by H.M. Smith 2018-04-20
    return (5*y+2)/(8*(1+y)**(5/2))

def Ka(y):
    #analytical form of H(y) - 3 F(y)/2
    return (2*y-1)/(8*(1+y)**(5/2))

def F2a(y):
    return 3*y/(8*(1+y)**(5/2))

def Ka2(y):
    #analytical form of H(y) - 3 F(y)/2
    return F2a(y) - Fa(y)/2

if __name__ == "__main__":
    # visualize the above functions
    # and test against known as asymptotes
    
    test1 = False
    test2 = False
    test3 = True
    import matplotlib.pyplot as plt
    import numpy as np

    me = 9.10938356e-31
    mC = 1.9926467051999998e-26
    mC = 1*mC
    y = me/mC;

    print("K, y>>1")
    print(K(1/y)*np.sqrt(mC) * (1 + 1/y))
    print(0.25*y**(3/2)*np.sqrt(mC) * (1 + 1/y))
    print("F, y>>1")
    print(F(1/y)*np.sqrt(mC) * (1 + 1/y))
    print(0.25*y**(3/2)*np.sqrt(mC) * (1 + 1/y))
    print("K, y<<1")
    print(K(y)*np.sqrt(me) * (1 + me/mC))
    print(-(1/8) *np.sqrt(me) * (1 + me/mC))
    print("F, y<<1")
    print(F(y)*np.sqrt(me) * (1 + me/mC))
    print(0.25*np.sqrt(me) * (1 + me/mC))

    if test1:
        y0 = np.linspace(0.5,4)
        y1 = np.linspace(0.001,0.01)
        y2 = np.linspace(100,1000)
        ones = np.ones(50)

        # since the functions are not properly
        # vectorized
        _F0 = [F(yi) for yi in y0]
        _H0 = [H(yi) for yi in y0]
        _K0 = [K(yi) for yi in y0]
        _F1 = [F(yi) for yi in y1]
        _H1 = [H(yi) for yi in y1]
        _K1 = [K(yi) for yi in y1]
        _F2 = [F(yi) for yi in y2]
        _H2 = [H(yi) for yi in y2]
        _K2 = [K(yi) for yi in y2]

        #p1=plt.plot(y,_F,'b',label='F')
        #p2=plt.plot(y,_H,'r',label='H')
        #p3=plt.plot(y,0.25/y,'b--',label=r'$\frac{1}{4y}$')
        #p4=plt.plot(y,0.25/y**2,'r--',label=r'$\frac{1}{4y^2}$')

        f, axarr = plt.subplots(3, sharex=True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel("y")

        axarr[0].plot(y0,_F0,'r',label='F')
        axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[0].legend()

        axarr[1].plot(y0,_H0,'r',label='H')
        axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[1].legend()

        axarr[2].plot(y0,_K0,'r',label='K')
        axarr[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[2].legend()

        plt.savefig("FHK.png")


        # 0 asymptote

        f, axarr = plt.subplots(3, sharex=True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel("y")

        axarr[0].plot(y1,_F1,'r',label='F')
        axarr[0].plot(y1,0.25*ones,'k--',label=r'$\frac{1}{4}$')
        axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[0].legend()

        axarr[1].plot(y1,_H1,'r',label='H')
        axarr[1].plot(y1,0.25*ones,'k--',label=r'$\frac{1}{4}$')
        axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[1].legend()

        axarr[2].plot(y1,_K1,'r',label='K')
        axarr[2].plot(y1,-(1/8)*ones,'k--',label=r'$-\frac{1}{8}$')
        axarr[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[2].legend()

        plt.savefig("FHKasym1")

        # infinity asymptote

        f, axarr = plt.subplots(3, sharex=True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel("y")

        axarr[0].plot(y2,_F2,'r',label='F')
        axarr[0].plot(y2,1/(4*y2**1.5),'k--',label=r'$\frac{1}{4y^{3/2}}$')
        axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[0].legend()

        axarr[1].plot(y2,_H2,'r',label='H')
        axarr[1].plot(y2,5/(8*y2**1.5),'k--',label=r'$\frac{5}{8y^{3/2}}$')
        axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[1].legend()

        axarr[2].plot(y2,_K2,'r',label='K')
        axarr[2].plot(y2,1/(4*y2**1.5),'k--',label=r'$\frac{1}{4y^{3/2}}$')
        axarr[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[2].legend()

        plt.savefig("FHKasym2")
        
    if test2:
        y0 = np.linspace(0.5,4)

        # since the functions are not properly
        # vectorized
        _F0 = [F(yi) for yi in y0]
        _H0 = [H(yi) for yi in y0]
        _K0 = [K(yi) for yi in y0]


        Fa0 = Fa(y0)
        Ha0 = Ha(y0)
        Ka0 = Ka(y0)

        
        f, axarr = plt.subplots(3, sharex=True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel("y")

        axarr[0].plot(y0,_F0,'r',label='F')
        axarr[0].plot(y0,Fa0,'k--',label='Fa')
        axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[0].legend()

        axarr[1].plot(y0,_H0,'r',label='H')
        axarr[1].plot(y0,Ha0,'k--',label=r'Ha')
        axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[1].legend()

        axarr[2].plot(y0,_K0,'r',label='K')
        axarr[2].plot(y0,Ka0,'k--',label=r'Ka')
        axarr[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[2].legend()

        plt.savefig("FHK_analytic.pdf")
    if test3:
        y0 = np.linspace(0.5,4)

        # since the functions are not properly
        # vectorized
        _FC1 = [FC1(yi) for yi in y0]
        _F1 = [F(yi) for yi in y0]
        _FC2 = [FC2(yi) for yi in y0]

        F1 = Fa(y0)
        F2 = F2a(y0)

        
        f, axarr = plt.subplots(2, sharex=True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlabel("y")

        axarr[0].plot(y0,_FC1,'r',label='FC1')
        axarr[0].plot(y0,F1,'k--',label='F1')
        axarr[0].plot(y0,_F1,'b:',label='F')
        axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[0].legend()

        axarr[1].plot(y0,_FC2,'r',label='FC2')
        axarr[1].plot(y0,F2,'k--',label=r'F2')
        axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[1].legend()

       
        plt.savefig("F1F2_Chandrasehkar.pdf")

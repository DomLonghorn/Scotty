import numpy as np
from scipy import constants as constants
from scipy import interpolate as interpolate
from scipy import optimize as optimize
from scipy import special as special
from scipy import integrate as integrate
import sys



#Cyclotron & Plasma frequency functions
def wce_func(B):
    return constants.e*B/constants.m_e
def wpe_func(n):
    return np.sqrt(n*(constants.e)**2./(constants.m_e*constants.epsilon_0))


def vte_half_func(Te):
    return np.sqrt(1000.*Te*constants.e/constants.m_e)

def Dtensor_func(kperp,kpar,w,K):
    #dispersion components in Stix frame
    npar = kpar*constants.c/w
    nperp = kperp*constants.c/w

    N2 = npar**2. + nperp**2
    Nx = nperp
    Ny = 0.
    Nz = npar
    N = np.array([Nx,Ny,Nz])
    D = np.outer(N,N) - N2*np.identity(3)
    return K + D

def ES_disp(kperp,kpar,w,B,ne,Te,lmax=50,symlog = False):
    #print(nperp)
    try:
        if np.shape(kperp)[0]==2:
            kperp = kperp[0]+1j*kperp[1]
    except:
        pass
    k2 = kperp**2 + kpar**2
    N2 = k2*((constants.c/w)**2)
    #k2 = kperp*np.conjugate(kperp) + kpar*np.conjugate(kpar)
    wpe = wpe_func(ne)
    wce = wce_func(B)
    vte = vte_half_func(Te)
    
    b = (kperp*vte/wce)**2

    
    #kpar = npar*w/cn.c
    D = 0. + 1j*0
    lterms = int(lmax)
    l_array = range(-lterms,lterms+1)
    
    for l,ll in enumerate(l_array):
        xi = get_xi(w,kpar,ll,wce,vte)/np.sqrt(2)
        Zn = get_Z(xi)
        In = special.iv(ll,b)
        D += In*Zn
    D = np.exp(-b)*D*(w/(np.sqrt(2)*abs(kpar)*vte))
        
    try:
        Det = N2*(1 + (((wpe/vte)**2)/k2) * (1 + D))
    except OverflowError:
        Det = float('inf')
        
    if symlog: #symlog sometimes preferable for rootfinder
        return symlog(Det)
    else:
        return Det

#function: calculate plasma dispersion (Fried-Conte) function
def get_Z(x):
    #plasma dispersion function zeta(x)
    #if x > 1e3: #large arg asymptote.
    #    return (-1/x)*(1 + 1/(2*(x**2)))
    #else:
    #    return 1j*np.sqrt(np.pi)*special.wofz(x)
    return 1j*np.sqrt(np.pi)*special.wofz(x)

#function: calculate first derivative of Fried-Conte function.
def get_Zprime(x):
    return -2. * (1 + x * get_Z(x))

#function: calculate "xi", as used in hot plasma dispersion
#inputs: wce = eB/me, and vte = sqrt(Te/me)
def get_xi(w,kpar,l,wce,vte):
    #ksi: argument in hot dispersion relation
    #return (w + l*wce)/(np.abs(kpar)*(vte**2))
    return (w - l*np.abs(wce))/(np.abs(kpar)*(vte))
import numpy as np
from scipy import constants as constants
from scipy import interpolate as interpolate
from scipy import optimize as optimize
from scipy import special as special
from scipy import integrate as integrate
import hotplasmaconstants as cn
import mpmath
import cmath
import sys


Te = 2.5 #in keV
#Cyclotron & Plasma frequency functions
def wce_func(B):
    return constants.e*B/constants.m_e
def wpe_func(n):
    return np.sqrt(n*(constants.e)**2./(constants.m_e*constants.epsilon_0))

def vte_func(Te):
    return np.sqrt(2000.*Te*constants.e/constants.m_e)
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

#________________EM non-relativistic dispersion_______________
#function: dispersion relation for hot non-rel Maxwellian plasma
#From Swanson.
#inputs: ne in m^-3, Te in keV. lmax is max harmonic terms to compute
#output: scalar complex Hamiltonian
def hot_EM_nonrel_disp(kperp,kpar,w,B,ne,Te,lmax=50,symlog = False, give_tensor = False):
    #using dielectric tensor from J. Decker PhD thesis Sec. 2.1.4
    wce = wce_func(B)
    wpe = wpe_func(ne)
    vte = vte_func(Te)
    
    K = get_Ktensor(kperp,kpar,vte,wpe,wce,w,lmax)
    D = Dtensor_func(kperp,kpar,w,K)
    
    try: #needed because root-finder can sometimes mess up array dimensions
        Dnew = np.zeros((3,3),dtype=complex)
        Dnew[0,0] = D[0][0]
        Dnew[0,1] = D[0][1]
        Dnew[0,2] = D[0][2]
        
        Dnew[1,0] = D[1][0]
        Dnew[1,1] = D[1][1]
        Dnew[1,2] = D[1][2]
        
        Dnew[2,0] = D[2][0]
        Dnew[2,1] = D[2][1]
        Dnew[2,2] = D[2][2]
        
        D = Dnew*1.
    except:
        pass
    
    #instead of calculating the Hamiltonian, return the dispersion tensor
    if give_tensor==True:
        return D
    
    Det = np.linalg.det(D)
    if symlog: #symlog sometimes preferable for rootfinder
        return symlog(Det)
    else:
        return Det

#function: Dielectric tensor components of hot non-rel Maxwellian plasma
#From Swanson.
#inputs: vte = sqrt(2Te/me), wce = eB/me, lmax is max harmonic terms to compute
#output: K tensor    
def get_Ktensor(kperp,kpar,vte,wpe,wce,w,lmax):
    #dielctric components of hot EM disp. relation for maxwellian plasma. 
    #Non-relativistic
    #From Laqua. Rewritten from Stix.
    mu = 0.5*((kperp*vte/wce)**2)
    
    l = 0
    T = get_Ttensor(kperp,kpar,vte,wce,mu,l,w)

    tol = 1e-20
    change = 10
    
    while change > tol and l < lmax:
        l +=1
        minT = get_Ttensor(kperp,kpar,vte,wce,mu,-l,w)
        plusT = get_Ttensor(kperp,kpar,vte,wce,mu,l,w)
        change = np.linalg.norm(minT+plusT)/np.linalg.norm(T)
        T += minT + plusT
        if l > lmax:
            change = tol/10.
    
    if mu.real < 700:
        K = np.identity(3,dtype=complex) + T * (wpe**2) * np.exp(-mu) / (w * kpar)
    else:
        return np.identity(3,dtype=complex)

    return K
# def get_Ktensor(kperp,kpar,vte,wpe,wce,w,lmax):
#     #dielctric components of hot EM disp. relation for maxwellian plasma. 
#     #Non-relativistic
#     #From Laqua. Rewritten from Stix.
#     mu = 0.5*((kperp*vte/wce)**2)
    
#     #kpar = npar*w/cn.c
    
#     #lterms = int(lmax)
#     #l_array = range(-lterms,lterms+1)
    
#     T = np.zeros((3,3),dtype=complex)
    
#     tol = 1e-15
#     change = 10
    
#     l = 0
#     while change > tol:
#         minT = get_Ttensor(kperp,kpar,vte,wce,mu,-l,w)
#         #print('shape',np.shape(minT))
#         plusT = get_Ttensor(kperp,kpar,vte,wce,mu,l,w)
#         change = np.linalg.norm(minT+plusT)/np.linalg.norm(T)
#         T += minT + plusT
#         l +=1
#         if l > lmax:
#             change = tol/10.
#     #print('last l',l)
    
#     if mu.real < 700:
#         K = np.identity(3,dtype=complex) + T * (wpe**2) * np.exp(-mu) / (w * kpar)
#     else:
#         return np.identity(3,dtype=complex)

#     #print('T',T)
#     #print('coef',xi0*(wpe_tmp/w)**2)
#     return K

#function: compute T tensor used in hot non-rel dielectric
#From Swanson
#inputs: vte = sqrt(2Te/me), mu = (kperp*vte/wce)**2, wce = eB/me, l is harmonic, w is 2*pi*f
#output: T tensor
def get_Ttensor(kperp,kpar,vte,wce,mu,l,w):
    try:
        xi = get_xi(w,kpar,l,wce,vte)*np.sign(kpar.real)
        #print('xi',xi)
        Zn = get_Z(xi)
        Znp = get_Zprime(xi)
        In = special.iv(l,mu)
        Inp = special.ivp(l,mu)
        sigma_e = +1
        
        T0 = 2 * mu * (In - Inp) * Zn / vte
        T1 = (l**2) * In * Zn / (mu * vte)
        T2 = 1j*sigma_e * l * (In - Inp) * Zn / vte
        T3 = - In * xi * Znp / vte
        T4 = kperp * l * In * Znp / (2*wce*mu)
        T5 = 1j * kperp * sigma_e * (In - Inp) * Znp / (2*wce)
        
        T = np.zeros((3,3),dtype=complex)
        T[0,0] = T1
        T[0,1] = T2
        T[0,2] = T4
        T[1,0] = -T[0,1]
        T[1,1] = T1 + T0
        T[1,2] = -T5
        T[2,0] = 1*T[0,2]
        T[2,1] = -T[1,2]
        T[2,2] = T3
        return T
    except OverflowError:
        return float('inf') 
    
def findX(wce, w):
    return 1 - (wce**2)/(w**2)

def Muller_method_rootfinder(guess,kpar,w,B,ne,Te,func,lmax = 50,roottol=1e-8,\
                             maxsteps=20,scl_guess=0.1, verbose = False):
    #"verbose = True" used for debugging
    istop = 0
    w_c = w/constants.c
    
    kguess = guess
    #findroot function requires three initial guesses
    #we take the provided guess and take neighboring points above and below
    #keep imaginary part fixed for the three guesses
    kguess_array = ((1-scl_guess)*kguess.real + 1j*kguess.imag, kguess, (1+scl_guess)*kguess.real + 1j*kguess.imag)
    
    try:
        kperp = complex(\
                    mpmath.findroot(
                        lambda kperp: func(float(kperp.real) + 1j*float(kperp.imag),kpar,w,B,ne,Te,lmax = lmax),\
                        kguess_array, solver='muller',tol=roottol,maxsteps=maxsteps, verbose = verbose)\
                        )
        f = func(kperp*w_c,kpar,w,B,ne,Te,lmax,False)
    except:
        kperp = np.nan
        f = np.nan
        istop = 1
        pass
    
    return kperp


def cold_disp_func(n,B,w,npar): #f isn't angular
    #cold dispersion relation
    #P,S,D are Stix dielectric components


    if B > 0:
        sign = 1.
    else:
        sign = -1.
        B = abs(B)
    wpe = np.sqrt(n*(cn.Ze*cn.e)**2./(cn.me*cn.eps))
    wpi = np.sqrt(n*(1.-cn.impfrac)*(cn.Zi*cn.e)**2/(cn.mi*cn.eps))
    wpim = np.sqrt(n*(cn.impfrac)*(cn.Zim*cn.e)**2/(cn.mim*cn.eps))
    wce = cn.Ze*cn.e*B/cn.me
    wci = cn.Zi*cn.e*B/cn.mi
    wcim = cn.Zim*cn.e*B/cn.mim
    #print('w','wpe','wce','wci',w,wpe,wce,wci)
    #print('ratio w/wLH', w/np.sqrt(abs(wce)*wci))
    #print('wpim,wcim',wpim,wcim)
    P =(1 - wpe**2./w**2
             - wpi**2/w**2
             - wpim**2/w**2)

    S = (1-wpe**2./(w**2-wce**2)-wpi**2/(w**2-wci**2)-
            wpim**2./(w**2-wcim**2))
    D = (wce*wpe**2./(w*(w**2-wce**2)) + wci*wpi**2/(w*(w**2-wci**2)) +
            wcim*wpim**2./(w*(w**2-wcim**2)))
    #print('approx D in LH range',abs(wpe**2/(w*wce)))
    #print('P,S,D',P,S,D)
    return P,S,sign*D

#function: calculate nperp from cold plasma dispersion
#input: dielectric terms in Stix notation, and npar
#return: nperp
def cold_nperp_func(P,S,D,npar,mode):
    #calc nperp = nperp(P,S,D) from Stix cold dispersion relation
    P0 = P*((npar**2.-S)**2-D**2)
    P2 = (S + P)*(npar**2. - S) + D**2
    P4 = S

    delta = (P2)**2 - 4.*P0*P4
    nperp2s = (-P2+cmath.sqrt(delta))/(2.*P4)
    nperp2f = (-P2-cmath.sqrt(delta))/(2.*P4)

    #The highest n_perp is always to slow mode

    if nperp2f > nperp2s:
        temp = nperp2f
        nperp2f = nperp2s
        nperp2s = temp
    if mode == -1:
        return cmath.sqrt(nperp2f)
    elif mode == 1:
        return cmath.sqrt(nperp2s) #SLOW O mode for current case SLOW x-mode above UHR
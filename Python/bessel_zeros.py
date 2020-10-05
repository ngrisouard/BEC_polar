# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:41:59 2020

@author: rundo
"""
import numpy as np
from scipy import special

def fi(y):
    c1 = 1.570796
    if y == 0:
        return 0
    elif y > 1e5:
        return c1
    else:
        if y<1:
            p = (3*y) ** (1/3)
            pp = p**2
            p = p * (1 + pp*(pp*(27 - 2*pp) - 210)/1575)
        else:
            p = 1/(y + c1)
            pp = p**2
            p = c1 - p*(1 + pp*(2310 + pp*(3003 + pp*(4818 + pp*(8591 + pp*16328))))/3465)
        pp = (y + p)**2
        r = (p - np.arctan(p+y)) / pp
        return p - (1+pp)*r*(1 + r/(p+y))
    
def bessr(d, a, x):
    if d == 1:
        return special.jv(a, x) / special.jv(a+1, x)
    elif d == 2:
        return special.yv(a, x) / special.yv(a+1, x)
    elif d == 3:
        return a/x - special.jv(a+1, x) / special.jv(a, x)
    else:
        return a/x - special.yv(a+1, x) / special.yv(a, x)
    
def bessel_zeros(d, a, n, e):
    
    z = np.zeros(n)
    
    aa = a**2
    mu = 4*aa
    mu2 = mu**2
    mu3 = mu**3
    mu4 = mu**4
    
    if d < 3:
        p = 7*mu - 31
        p0 = mu - 1
        if (1+p) == p:
            p1 = 0
            q1 = 0
        else:
            p1 = 4*(253*mu2 - 3722*mu+17869)*p0/(15*p)
            q1 = 1.6*(83*mu2 - 982*mu + 3779)/p
    
    else:    
        p = 7*mu2 + 82*mu - 9
        p0 = mu + 3
        if (p+1) == 1:
            p1 = 0
            q1 = 0
        else:
            p1 = (4048*mu4 + 131264*mu3 - 221984*mu2 - 417600*mu + 1012176)/(60*p)
            q1 = 1.6*(83*mu3 + 2075*mu2 - 3039*mu + 3537)/p
            
    if d == 1 or d == 4:
        t = 0.25
    else:
        t = 0.75
    tt = 4*t
    
    if d < 3:
        pp1 = 5/48
        qq1 = -5/36
        
    else:
        pp1 = -7/48
        qq1 = 35/288
    
    y = 0.375 * np.pi
    
    if a >= 3:
        bb = a ** (-2/3)
    else:
        bb = 1
    a1 = 3*a - 8
    
    for s in range(1, n+1):
        if a == 0 and s == 1 and d == 3:
            x = 0
            j = 0
        else:
            if s >= a1:
                b = (s + 0.5*a - t) * np.pi
                c = .015625 / (b**2)
                x = b - 0.125*(p0 - p1*c)/(b*(1 - q1*c))
            else:
                if s == 1:
                    if d == 1:
                        x = -2.33811
                    elif d == 2:
                        x = -1.17371
                    elif d == 3:
                        x = -1.01879
                    else:
                        x = -2.29444
                else:
                    x = y*(4*s - tt)
                    v = x ** (-2)
                    x = -x**(2/3) * (1 + v*(pp1 + qq1*v))
                
                u = x*bb
                v = fi(2/3 * (-u)**1.5)
                w = 1/np.cos(v)
                xx = 1 - w**2
                c = np.sqrt(u/xx)
                
                if d < 3:
                    x = w*(a + c*(-5/u - c*(6 - 10/xx))/(48*a*u))
                else:
                    x = w*(a + c*(7/u + c*(18 - 14/xx))/(48*a*u))
            j = 0
        
            while ((j==0) or ((j<5) and (abs(w/x)>e))):
                xx = x**2
                x4 = x**4
                a2 = aa - xx
                r0 = bessr(d, a, x)
                j = j+1
                if d < 3:
                    u = r0
                    w = 6*x*(2*a + 1)
                    p = (1 - 4*a2)/w
                    q = (4*(xx-mu) - 2 - 12*a)/w
                else:
                    u = -xx*r0/a2
                    v = 2*x*a2/(3*(aa+xx))
                    w = 64*a2**3
                    q = 2*v*(1 + mu2 + 32*mu*xx + 48*x4)/w
                    p = v*(1 + (40*mu*xx + 48*x4 - mu2)/w)
                
                w = u*(1 + p*r0)/(1 + q*r0)
                x = x+w
        
            z[s-1] = x
    return z
    
    
        
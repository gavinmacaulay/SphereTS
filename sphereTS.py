# -*- coding: utf-8 -*-
"""
Created on Thu May 01 14:12:05 2014

@author: gavinj
"""
from __future__ import division
import math
import cmath
from scipy import special

from numpy import arange

def sphereTS(f, a, c, c1, c2, rho, rho1):
    """
    Calculates the acoustic backscattered target strength of an elastic sphere immersed in a fluid.
    
    Implements equations 6 to 9 in MacLennan (1981). 
    
    Args:
        :f: acoustic frequency [Hz]
        :a: sphere radius [m]
        :c: sound speed in the fluid that surrounds the sphere [m/s]
        :c1: longitudinal sound speed of the sphere material [m/s]
        :c2: transverse sound speed of the sphere material [m/s]
        :rho: density of the fluid that surrounds the sphere [kg/m^3]
        :rho1: density of the sphere material [kg/m^3]

    Returns:
        The target strength of the sphere [dB re 1m^2]
      
    Raises:
        ValueError if any of the input arguments are less than zero.        
    
    References:
        MacLennan, D.N., 1981. The Theory of Solid Spheres 
        as Sonar Calibration Targets. Scottish Fisheries Research 
        Report No. 22, Department of Agriculture and Fisheries for Scotland. 
    """
    
    if f < 0:  
        raise ValueError("Frequency must be greater than 0 (was {}).".format(f))
    if a < 0:
        raise ValueError('Sphere radius must be greater than 0 (was {}).'.format(a))
    if c < 0:
        raise ValueError('Fluid sound speed must be greater than 0 (was {})'.format(c))
    if c1 < 0:
        raise ValueError('Longitudinal sound speed of the sphere must be greater than 0 (was {})'.format(c1))
    if c2 < 0:
        raise ValueError('Transverse sound speed of the sphere must be greater than 0 (was {})'.format(c2))
    if rho < 0:
        raise ValueError('Density of the fluid must be greater than 0 (was {})'.format(rho))
    if rho1 < 0:
        raise ValueError('Density of the sphere must be greater than 0 (was {})'.format(rho1))
    
    k = 2.0*math.pi*f/c
    q = k*a
    q1 = q*c/c1
    q2 = q*c/c2
    alpha = 2.0 * (rho1/rho) * ((c2/c)**2)
    beta = (rho1/rho) * ((c1/c)**2) - alpha 
    form = 0.0
    tol = 1e-10
    pi2 = math.pi / 2.0
    
    for l in range(0, 50):
        # The relationship for the derivatives of the bessel
        # functions are from Rudgers, A.J., 1967. Techniques 
        # for Numerically Evaluating the Formulas Describing 
        # Monostatic Reflection of Acoustic Waves by Elastic 
        # Spheres. NRL Report 6551, Acoustic Research Branch, 
        # Sound Division, Naval Research Laboratory.
        j_q =    special.jv(l+0.5, q) * math.sqrt(pi2/q)
        jm1_q =  special.jv(l-0.5, q) * math.sqrt(pi2/q)
        j_qd =   jm1_q - (l+1) / q * j_q
        j_q1 =   special.jv(l+0.5, q1) * math.sqrt(pi2/q1)
        jm1_q1 = special.jv(l-0.5, q1) * math.sqrt(pi2/q1)
        j_q1d =  jm1_q1 - (l+1) / q1 * j_q1
        j_q1dd = 1 / (q1**2) * ((l+1)*(l+2) - q1**2) * j_q1 - 2.0 / q1 * jm1_q1
        j_q2 =   special.jv(l+0.5, q2) * math.sqrt(pi2/q2)
        jm1_q2 = special.jv(l-0.5, q2) * math.sqrt(pi2/q2)
        j_q2d =  jm1_q2 - (l+1) / q2 * j_q2
        j_q2dd = 1 / (q2**2) * ((l+1)*(l+2) - q2**2) * j_q2 - 2.0 / q2 * jm1_q2
       
        y_q =   special.yv(l+0.5, q) * math.sqrt(pi2/q)
        yp1_q = special.yv(l+1.5, q) * math.sqrt(pi2/q)
        y_qd = l / q * y_q - yp1_q
       
        A2 = (l**2 + l-2.0)*j_q2 + (q2**2)*j_q2dd
        A1 = 2.0*l*(l+1) * (q1 * j_q1d - j_q1)
        B2 = A2 * q1**2 * (beta*j_q1-alpha*j_q1dd) - A1*alpha*(j_q2-q2*j_q2d)
        B1 = q * (A2 * q1 * j_q1d - A1 * j_q2)
        neta = math.atan(-(B2*j_qd - B1*j_q)/(B2*y_qd - B1*y_q))
    
        newterm = (-1.0)**l * (2.0*l+1) * cmath.sin(neta) * cmath.exp(1j * neta)
        
        form = form + newterm
        if abs(newterm)/abs(form) < tol:
            break
    
    form = -2.0/q*form
    sigma = math.pi * a**2 * abs(form)**2
    return 10.0 * math.log10(sigma / (4.0*math.pi))

def sphereTSFreqResponse():
    """
    
    """
    a = 0.0381/2
    c = 1490
    c1 = 6853
    c2 = 4171
    rho = 1026.2
    rho1 = 14900
    fstart = 12e3
    fstop = 200e3
    
    for f in arange(fstart, fstop, 100):
        ts = sphereTS(f, a, c, c1, c2, rho, rho1)
        print('TS = {} at {} kHz.'.format(ts, f/1000.))
    pass

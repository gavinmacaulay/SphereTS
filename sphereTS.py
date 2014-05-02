# -*- coding: utf-8 -*-
"""
Created on Thu May 01 14:12:05 2014

@author: gavinj
"""
from __future__ import division
import math
import cmath
from scipy import special

def sphereTS(f, a, c, c1, c2, rho, rho1):
    """
    Calculates the backscattered acoustic TS of an elastic sphere in a fluid.
    
    Calculates the TS as defined by equations 6 to 9 in: 
    MacLennan, D.N., 1981. The Theory of Solid Spheres 
    as Sonar Calibration Targets. 22, Department of 
    Agriculture and Fisheries for Scotland. 
    
    Input variables are:
    f - acoustic frequency, [Hz]
    a - sphere radius [m]
    c - sound speed in water [m/s]
    c1 - longitudinal sound speed of the sphere material [m/s]
    c2 - transverse sound speed of the sphere material [m/s]
    rho - density of fluid that surrounds the sphere [kg/m^3]
    rho1 - density of sphere material [kg/m^3]
     
    The input variables should all be scalars.
    """
    
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

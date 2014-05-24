# -*- coding: utf-8 -*-
"""    
    This module provides a function that calculates the target strength of an 
    elastic sphere immersed in a fluid. 
    
    Addition functions are included to calculate the target strength at a 
    range of frequencies, provide default material properties for selected
    elastic materials, and to calculate seawater properties from measurements
    of temperature, salinity, and depth.    


    Copyright 2014 Gavin Macaulay

    This file is part of SphereTS.

    SphereTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SphereTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SphereTS.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division

import math
import cmath

from scipy import special
import numpy as np
import gsw

def sphereTS(f, a, c, c1, c2, rho, rho1):
    """
    Calculates the acoustic backscattered target strength of an 
    elastic sphere immersed in a fluid.
    
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
    
    if f <= 0:  
        raise ValueError("Frequency must be greater than 0 (was {}).".format(f))
    if a <= 0:
        raise ValueError('Sphere radius must be greater than 0 (was {}).'.format(a))
    if c <= 0:
        raise ValueError('Fluid sound speed must be greater than 0 (was {})'.format(c))
    if c1 <= 0:
        raise ValueError('Longitudinal sound speed of the sphere must be greater than 0 (was {})'.format(c1))
    if c2 <= 0:
        raise ValueError('Transverse sound speed of the sphere must be greater than 0 (was {})'.format(c2))
    if rho <= 0:
        raise ValueError('Density of the fluid must be greater than 0 (was {})'.format(rho))
    if rho1 <= 0:
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
    # precalculate some things
    sqrpi2q = math.sqrt(pi2/q)
    sqrpi2q1 = math.sqrt(pi2/q1)
    sqrpi2q2 = math.sqrt(pi2/q2)
    
    for l in range(0, 50):
        # The relationship for the derivatives of the bessel
        # functions are from Rudgers, A.J., 1967. Techniques 
        # for Numerically Evaluating the Formulas Describing 
        # Monostatic Reflection of Acoustic Waves by Elastic 
        # Spheres. NRL Report 6551, Acoustic Research Branch, 
        # Sound Division, Naval Research Laboratory.
        j_q =    special.jv(l+0.5, q) * sqrpi2q
        jm1_q =  special.jv(l-0.5, q) * sqrpi2q
        j_qd =   jm1_q - (l+1) / q * j_q
        j_q1 =   special.jv(l+0.5, q1) * sqrpi2q1
        jm1_q1 = special.jv(l-0.5, q1) * sqrpi2q1
        j_q1d =  jm1_q1 - (l+1) / q1 * j_q1
        j_q1dd = 1 / (q1**2) * ((l+1)*(l+2) - q1**2) * j_q1 - 2.0 / q1 * jm1_q1
        j_q2 =   special.jv(l+0.5, q2) * sqrpi2q2
        jm1_q2 = special.jv(l-0.5, q2) * sqrpi2q2
        j_q2d =  jm1_q2 - (l+1) / q2 * j_q2
        j_q2dd = 1 / (q2**2) * ((l+1)*(l+2) - q2**2) * j_q2 - 2.0 / q2 * jm1_q2
       
        y_q =   special.yv(l+0.5, q) * sqrpi2q
        yp1_q = special.yv(l+1.5, q) * sqrpi2q
        y_qd = l / q * y_q - yp1_q
       
        A2 = (l**2 + l-2.0)*j_q2 + (q2**2)*j_q2dd
        A1 = 2.0*l*(l+1) * (q1 * j_q1d - j_q1)
        B2 = A2 * q1**2 * (beta*j_q1-alpha*j_q1dd) - A1*alpha*(j_q2-q2*j_q2d)
        B1 = q * (A2 * q1 * j_q1d - A1 * j_q2)
        neta = math.atan(-(B2*j_qd - B1*j_q)/(B2*y_qd - B1*y_q))
    
        newterm = (-1.0)**l * (2.0*l+1) * cmath.sin(neta) * cmath.exp(1j * neta)
        
        form = form + newterm
        if abs(newterm/form) < tol:
            break
    
    form = -2.0/q*form
    sigma = math.pi * a**2 * abs(form)**2
    
    return 10.0 * math.log10(sigma / (4.0*math.pi))

def sphereTSFreqResponse(fstart, fstop, a, c, c1, c2, rho, rho1, fstep=100):
    """
    Calculates the acoustic backscattered target strength of an 
    elastic sphere immersed in a fluid at a range of frequencies.
    
    Calls sphereTS() to calculate the TS at discrete frequencies.
    
    Args:
        :fstart: lowest acoustic frequency [Hz]
        :fstop: highest acoustic frequency
        :fstep: frequency step size. Optional, defaults to 100 [Hz]
        :a: sphere radius [m]
        :c: sound speed in the fluid that surrounds the sphere [m/s]
        :c1: longitudinal sound speed of the sphere material [m/s]
        :c2: transverse sound speed of the sphere material [m/s]
        :rho: density of the fluid that surrounds the sphere [kg/m^3]
        :rho1: density of the sphere material [kg/m^3]
        
    Returns:
        :f: a NumPy array containing the frequencies at which the TS has 
            been calculated [Hz]
        :TS: a NumPy array containing the sphere TS at the frequencies 
             given in f [dB re 1 m^2]
    """
    
    f = np.arange(fstart, fstop, fstep)
    TS = np.zeros(len(f))
    
    for i in range(len(f)):
        TS[i] = sphereTS(f[i], a, c, c1, c2, rho, rho1)
        
    return f,TS
    
def materialProperties():
    """
    Provides sound speed and density values for selected sphere materials.
    
    Returns:
        A dictionary with keys containing the material name, and values being 
        another dictionary with keys 'rho1', 'c1', and 'c2' corresponding to
        density, transveral sound speed, and longitudinal sound speed 
        respectively. Units are kg/m^3, m/s, and m/s respectively.
    """
    
    return {'Tungsten carbide': {'c1': 6853.0, 'c2': 4171.0, 'rho1': 14900.0},
            'Copper':           {'c1': 4760.0, 'c2': 2288.5, 'rho1':  8947.0}, 
            'Stainless steel':  {'c1': 5610.0, 'c2': 3120.0, 'rho1':  7800.0},
            'Alumnium':         {'c1': 6260.0, 'c2': 3080.0, 'rho1':  2700.0}}

def calcWaterProperties(SP, t, p):
    """
    Calculates seawater density and sound speed using the TEOS-10 equations.
    
    Uses the gsw toolbox, which is a wrapper around the C-based TEOS-10 library.
    
    Args:
        :SP: Practical salinity [PSU]
        :t: Temperature [degC]
        :p: Pressure [dbar]
        
    Returns:
        :c: speed of sound in water [m/s]
        :rho: density of water [kg/m^3]
    """
    
    SR = gsw.SR_from_SP(SP)
    CT = gsw.CT_from_t(SR, t, p)
    
    c = gsw.sound_speed(SR, CT, p)
    rho = gsw.rho(SR, CT, p)
    
    return c, rho


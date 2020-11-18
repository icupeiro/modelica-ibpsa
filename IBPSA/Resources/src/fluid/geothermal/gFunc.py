# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:14:43 2018

@author: u0112721
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import datetime
import pygfunction as gt
from scipy.constants import pi
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import j0, j1, y0, y1, exp1
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd


def shortTermCorrection(time, gFunc, r_b, aSoi):
    
    def _CHS(u, Fo, p):
        CHS_integrand = 1./(u**2*pi**2)*(np.exp(-u**2*Fo) - 1.0) / (j1(u)**2 + y1(u)**2) * (j0(p*u)*y1(u) - j1(u)*y0(p*2))
        return CHS_integrand
    
    def _ILS(t, aSoi, dis):
        ILS = exp1(dis**2/(4*aSoi*t))
        return ILS

    for i in range(len(time)):
        ILS = _ILS(time[i], aSoi, r_b)
        CHS, err = quad(
            _CHS, 1e-12, 100., args=(aSoi*time[i]/r_b**2, 1.))
        gFunc[i] = gFunc[i] + 2*pi*CHS - 0.5*ILS

    return gFunc


def main():
    # -------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------

    # Borehole dimensions
    D = 10.0             # Borehole buried depth (m)
    H = [50, 100]           # Borehole length (m)
    r_b = 0.075         # Borehole radius (m)
    B = 6             # Borehole spacing (m)

    plt.rc('figure')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # Axis labels
    ax1.set_xlabel(r'$ln(t)$')
    ax1.set_ylabel(r'$g(t)/(2pi k_{soi} H)$')
    # Axis limits
    #ax1.set_xlim([-10.0, 5.0])
    #ax1.set_ylim([0., 20.])
    # Show minor ticks
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    # Adjust to plot window
    plt.tight_layout()
    lin = ['-', '-.']
    wid = [1, 2]
    i = 0
    j = -1

    for H in H:

        j = j + 1


        # Soil thermal properties
        kSoi = 1.5				# Ground thermal conductivity (W/(mK))
        cSoi = 1200				# Specific heat capacity of the soil (J/(kgK))
        dSoi = 1800				# Density of the soil (kg/m3)
        aSoi = kSoi/cSoi/dSoi   # Ground thermal diffusivity (m2/s)

        # Number of segments per borehole
        nSegments = 12

        # Geometrically expanding time vector. 

        nt = 100						   # Number of time steps
        ts = H**2/(9.*aSoi)            # Bore field characteristic time
        #ttsMax = np.exp(5)
        #ydes = 50						# Design projected period (years)
        dt = 3600.		                # (Control) Time step (s)
        tmax = ts*np.exp(5)                 # Maximum time
        #tmax = ttsMax*ts                # Maximum time

        time = gt.utilities.time_geometric(dt, tmax, nt)
        # -------------------------------------------------------------------------
        # Borehole fields
        # -------------------------------------------------------------------------

        # Field definition
        N_1 = 1
        N_2 = 1
        nBor = 1
        boreField = gt.boreholes.rectangle_field(N_1, N_2, B, B, H, D, r_b)

        gFunc = gt.gfunction.uniform_temperature(boreField, time, aSoi, nSegments=nSegments, disp=True)
        #gFunc = shortTermCorrection(time, gFunc, r_b, aSoi)
        gFunc = gFunc / (2*np.pi*kSoi*H*nBor) 

        kSoi2 = kSoi*2
        aSoi2 = kSoi2/cSoi/dSoi
        ts2 = H**2/(9.*aSoi2) 
        tmax2 = ts2*np.exp(5)                 # Maximum time


        time2 = gt.utilities.time_geometric(dt, tmax2, nt)    

        gFunc2 = gt.gfunction.uniform_temperature(boreField, time, aSoi2, nSegments=nSegments, disp=True)
        #gFunc2 = shortTermCorrection(time, gFunc2, r_b, aSoi2)
        gFunc2 = gFunc2 / (2*np.pi*kSoi2*H*nBor)


        #Adding zero as the first element
        time = np.insert(time, 0, 0)
        time2 = np.insert(time2, 0, 0)
        gFunc = np.insert(gFunc, 0, 0) 
        gFunc2 = np.insert(gFunc2,0,0)    

        # -------------------------------------------------------------------------
        # Figure
        # -------------------------------------------------------------------------


        # Draw g-function
        ax1.plot(np.log(time[1:]), gFunc[1:], color='blue', linestyle=lin[i], lw=1.5, label='k=1.5; H='+str(nBor*H))
        ax1.plot(np.log(time[1:]), gFunc2[1:], color='red', linestyle=lin[i], lw=1.5, label='k=3;H='+str(nBor*H))
        ax1.legend(loc='upper left')
        #ax1.axvline(x=np.log(6*30*24*60*60), color='black', linewidth=2)
        #ax1.axvline(x=np.log(3*30*24*60*60), color='black', linewidth=1)
        #ax1.axvline(x=np.log(6*30*24*60*60), color='red')
        #ax1.axvline(x=np.log(6*30*24*60*60), color='blue')
    
        i = 1
    plt.show()

    return



# Main function
if __name__ == '__main__':
    main()
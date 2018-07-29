# -*- coding: utf-8 -*-

import math

import numpy as np
from scipy import stats, signal


def uniform(x):
    """Uniform(square) kernel."""
    z = np.ones(len(x),np.float)
    z[np.abs(x) >= 1] = 0
    return z         
    
def triangle(x):
    """Triangular kernel."""
    z = (1.0 - np.abs(x))
    z[np.abs(x) >= 1] = 0
    return z        

def tricube(x):
    """Tricube kernel."""
    z = (1.0 - np.abs(x)**3)**3
    z[np.abs(x) >= 1] = 0
    return z
    
def cubed(x):
    """Cubic kernel."""
    z = (1.0 - np.abs(x)**3)
    z[np.abs(x) >= 1] = 0
    return z        

def triweight(x):
    """Triweight kernel."""
    z = 1.09375 * ((1.0 - np.array(x)**2)**3)
    z[np.abs(x) >= 1] = 0
    return z
    
def epanechnikov(x):
    """Epanechnikov kernel."""
    z = 0.75 * (1.0 - np.array(x)**2)
    z[np.abs(x) >= 1] = 0
    return z   
        
def quartic(x):
    """Quartic kernel."""
    z = 0.9375 * ((1.0 - np.array(x)**2)**2)
    z[np.abs(x) >= 1] = 0
    return z  
    
biweight = quartic 


rsqrt2pi = 1./math.sqrt(2.0*math.pi)

def gaussian(x):
    """Gaussian kernel."""
    z = rsqrt2pi * np.exp(-0.5 * x**2)
    return z    
    
    
pi4 = math.pi/4.
pi2 = math.pi/2.

def cosine(x):
    """Cosine kernel."""
    z = pi4 * np.cos(pi2*x)
    z[np.abs(x) >= 1] = 0
    return z    
        


def kernel_smooth(x,y,h, kernel=uniform, ctinband=False):
    """Calculcate a `kernel smoothed' moving average of y, at the given
    x-coords and with half-bandwith, h.
    
    This is a function to calculate moving averages given uneven sampling. 
    Calculates moving average of y at coordinate x_i by weighted averaging 
    over all points in range (x_i-h, x_i+h). Weights are given by the kernel
    function.
    
    This is equivalent to the Nadaraya-Watson regression estimate.
    
    """
    x = np.array(x)
    y = np.array(y)
    
    olen = len(x)
    
    # at the head and tail of the sequence reflect the right and left (respectively)
    # half-bandwidths
    xfirst = x[0]
    startband = x <= xfirst + h
    xstart = x[startband]
    ystart = y[startband]
    nstart = len(xstart)-1
    
    xlast = x[-1]
    endband = x >= xlast - h
    xend = x[endband]
    yend = y[endband]
    nend = len(xend) - 1    
    x = np.hstack((xstart[::-1][:-1], x, xend[::-1][1:]))
    y = np.hstack((ystart[::-1][:-1], y, yend[::-1][1:]))
    
    lx, ly = len(x), len(y)
    if lx != ly:
        raise Exception("x and y must be same length.")
    z = []
    ninband = []
    for i in range(lx):
        c = x[i]
        inband =  np.logical_and(x >= c-h, x <= c+h)
        if ctinband:
            ninband.append(len(np.flatnonzero(inband)))
        xfrac = (np.abs(x[inband] - c))/float(h)
        xwt = kernel(xfrac)
        ywin = np.sum(y[inband]*xwt)/np.sum(xwt)
        z.append(ywin)
    
    # trim off the reflected right/left half-bandwidths
    if ctinband:
        return z[nstart:olen+nstart], ninband[nstart:olen+nstart]        
    return z[nstart:olen+nstart]   





def connected_intervals(x, y, maxgap=10000):
    """Determine the connected intervals over a set of x,y observations by 
    identifying the 'gaps'.
    
    Gaps are defined as adjacent x values where x[i+1]-x[i] > maxgap
    
    This is useful when you want to draw a plot over a set of data but you 
    don't want to connect points that span intervals where there is no data.
    """
    x = np.array(x)
    y = np.array(y)
    diff = x[1:] - x[:-1]
    gaps = [i+1 for i in np.flatnonzero(diff > maxgap)]
    if not len(gaps):
        return [x],[y]
    newx, newy = [], []
    idx = 0
    for i,j in enumerate(gaps):
        jx = x[idx:j]
        jy = y[idx:j]
        newx.append(jx)
        newy.append(jy)
        if i == len(gaps)-1:
            newx.append(x[j:])
            newy.append(y[j:])
        idx = j
    return newx, newy
            


def minfilter(x, size=5):
    return signal.order_filter(x, np.ones(size), 0)

def maxfilter(x, size=5):
    return signal.order_filter(x, np.ones(size), size-1)
    #return signal.order_filter(x, np.ones(size), size-1)

def minmaxfilter(x, size=5):
    return minfilter(maxfilter(x, size), size)

def maxminfilter(x, size=5):
    return maxfilter(minfilter(x, size), size)

def connected_maxminfilter(x, y, size=5, maxgap=10000):
    gapx, gapy = connected_intervals(x, y, maxgap)
    maxminy = [maxminfilter(i, size) for i in gapy]
    return np.concatenate(gapx), np.concatenate(maxminy)

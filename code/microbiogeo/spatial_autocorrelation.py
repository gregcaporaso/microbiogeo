#!/usr/bin/env python
# File created on 8 Aug 2011

from __future__ import division
from numpy import tri, argsort, array, ones, arange, random, mean, unique, zeros, exp
from numpy.random import permutation
import sys
import os.path
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat, parse_mapping_file
from cogent.maths.stats.test import zprob

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"


def morans(w, x):
    """Computes Moran's I
    
       PARAMETERS
       w: a matrix of spatial weights
       x: a vector of variables of interested, taken from mapping file
       
       RETURNS
       I: value computed by Moran's I
       
    """
    # Average of x, used for moran's computations
    xbar = mean(x)
    
    # left side of the equation
    accumulator = 0
    for i in range(len(w)):
        for j in range(len(w)):
            accumulator = accumulator + w[i][j]
    left = len(x)/accumulator
    
    # right side of the equation
    accumulator = 0
    for i in range(len(x)):
        for j in range(len(x)):
            accumulator = accumulator + (w[i][j] * (x[i]-xbar) * (x[j]-xbar))

    # Top referes to the top half of the equation to calculate I
    # Bottom refers to the bottom half, 
    top = accumulator
    accumulator = 0    
    for i in range(len(x)):
        accumulator = accumulator + ((x[i] - xbar)**2)
    bottom = accumulator
    right = top/bottom

    I = left * right
    return I


def morans_variance(w, x):
    """Compute Moran's I and find varience
    
        PARAMETERS
        w: a matrix of spatial weights
        x: a vector of variables of interested, taken from mapping file
        
        RETURNS
        i_value: result of Moran's I
        v_value: calculated varience
    
    """
    # Calculate the I value
    i_value = morans(w, x)

    # Equation variables
    xbar = mean(x)
    n = len(x)
    sum_w = sum(sum(w))
    
    s1 = _s1(w)
    s2 = _s2(w)
    s3 = _s3(w, x, xbar, n)
    s4 = _s4(n, s1, s2, sum_w)
    s5 = _s5(s1, n, sum_w) 
    
    # Var(I)
    """
    The equations that should be used in this section to find the right
    values for varience were not resolved completely, the ones we have tried
    and suspect should work come from a book: "Spatial Processes Models &
    Applications" by AD Cliff & JK Ord. QA 278.2 C62
    ISBN 0 85086 081 4. Around pate 46
    
    The major question that needs to be answered is whether the equation should
    be the ones for assumption R or assumption N
    
    Equations for s1-s5 are taken from the wiki page for Moran's I.
    """
    top = (n * s4) - (s3 * s5)
    bottom = (n - 1) * (n - 2) * (n - 3) * (sum_w**2)
    v_value = top/bottom
    e_value = -1/(n-1)
    z_value = (i_value - e_value)/ v_value**(.5)
    e2 = (((n**2)*s1) - (n*s2) + (3*(sum_w**2)))/((n-1)*(n+1)*(sum_w**2))
#    m2 = 
#    e2r = (n(((n**2) - (3n) + 3)*s1)) - (n*s2) + (3*(sum_w**2))) - 
    print "e[i2]: ",e2
    print "e[i**2]-e[i]**2: ", (v_value - (e_value**2))
    print "e[i]**2: ",(e_value**2)
    print "i: ",i_value
    print "e: ",e_value
    print "v: ",v_value
    print "z: ",z_value
    print "p: ",zprob(z_value)
    return i_value, v_value


def _format_morans_results(input_path, i_value, v_value='NA'):
    """Format for output file
        
        First Line(header): Input_filepath  Morans I value  Varience
        Second Line: Values for the above headers
        Everything is tab delimited
     
    """
    result = ['Input_filepath\tMorans I value\tVariance']
    result.append('{0}\t{1}\t{2}'.format(input_path,i_value,v_value))
    return result


def _w_inverse(distmtx):
    """Weight w matrix using inverse method
        
        PARAMETERS
        distmtx: a distance matrix
        
        RETURNS:
        w: the distance matrix weighted by their inverses
    
    """
    # Set Diagonal to 1
    for i in range(len(distmtx)):
        distmtx[i][i] = 1
    
    # Inverse the matrix
    w = 1/distmtx # does not work in the case that there is a 0 value

    # Set Diagonal to 0
    for i in range(len(distmtx)):
        w[i][i] = 0
    
    return w
    
    
def _w_exponential(distmtx, kvalue):
    """Weight w matrix using exponential method
    
        PARAMETERS
        distmtx: a distance matrix
        kvalue: integer to use as the k value in the equation
        
        equation: e**(-k * w)

        RETURNS
        w: the distance matrix weighted by the exponential equation
    
    """
    # Calculate the values of w based on the given exp equation
    w = exp(-kvalue * distmtx)
    
    # Set diagonal to 0
    for i in range(len(distmtx)):
        w[i][i] = 0
    return w
  

# _s1 through _s5 are taken from the Wiki page for Moran's I

def _s1(w):
    """S1, part of the varience computation
    
        S1 = .5 * (SUMi(SUMj((w[i][j] + w[j][i])**2)))
    
    """
    accumulator = 0
    for i in range(len(w)):
        for j in range(len(w)):
            accumulator = accumulator + ((w[i][j] + w[j][i])**2)
            print "accumulator: ",accumulator
    return .5 * accumulator
    
    
def _s2(w):
    """S2, part of the varience computation
    
        S2 = SUMi(SUMj(w[i][j] + SUMj(w[j][i])))**2
    
    """
    accumulator1 = 0
    accumulator2 = 0
    accumulator3 = 0
    for i in range(len(w)):
        for j in range(len(w)):
            accumulator2 = accumulator2 + w[i][j]
        for j in range(len(w)):
            accumulator3 = accumulator3 + w[j][i]               
        accumulator1 = accumulator1 + ((accumulator2 + accumulator3)**2)
        accumulator2 = 0
        accumulator3 = 0
    return accumulator1
    
def _s3(w, x, xbar, n):
    """S3, part of the varience computation
    
        S3 = N**-1 * SUMi((x[i] - xbar)**4)
            (N**-1 * SUMi((x[i] - xbar)**2))**2
            
    """

    accumulator = 0
    for i in range(len(x)):
        accumulator = accumulator + ((x[i] - xbar)**4)
    top = (n**-1) * accumulator
    
    accumulator = 0
    for i in range(len(x)):
        accumulator = accumulator + ((x[i] - xbar)**2)
    bottom = ((n**-1) * accumulator)**2
    
    return top/bottom
    
    
def _s4(n, s1, s2, sum_w):
    """S4, part of the varience computation
    
        S4 = (N**2 - 3N + 3) * S1 - N * S2 + 3(SUMi(SUMj(w[i][j])))**2
    
    """
    return (((n**2 - (3 * n) + 3) * s1) - (n * s2) + (3*(sum_w**2)))/1
    
    
def _s5(s1, n, sum_w):
    """S5, part of the varience computation
    
        S5 = S1 - 2 * N * S1 + 6(SUMi(SUMj(w[i][j])))**2
        
    """
    return s1 - (2 * n * s1) + ((6 * (sum_w**2)))
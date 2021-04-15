# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 19:49:02 2021

@author: franc
"""
import numpy as np
import scipy.special as sp

#
def simpson(func, a, b, N):
    """
    Integration by Simpson method
    param 1 func: func is the function to integrate 
    param 2 a: a is the lower boundary of the integral 
    param 3 b: b is the upper boundary of the integral  
    param 4 N: N is the number of slic
    return: I, the integration of the function   
    """
    h = abs(b - a) / N  # Size of the slices 
    x = np.linspace(a, b, N+1) #N+1 points in the interval [a,b]
    y = func(x)
    
    y_odd = y[1:-1:2] #Sum of odd 
    s1 = np.sum(y_odd)
    
    y_even = y[2:-2:2] #Sum of even
    s2 = np.sum(y_even)
    
    I = h/3*(y[0] + y[-1] + s1*4+s2*2) #Area under the curve 
    
    return I 
    
def Hankel_transform(func, p, n, L):
    """
    Parameters
    ----------
    g : function that receives one argument and returns one 
        g the function to transform
    p :The p value in the Fourrier-Bessel space 
    n : integer
        The order of the Bessel function of the first kind.
    L : The number of term for the partial sum 
    Returns
    -------
    The partial sum L of the hankel transform of the function func 
    """
    zeros = sp.jn_zeros(n, L+1) # array of the L+1 first zeros of bessel order n  
    value = 0 #Initialize the value 
    
    #function to integrate
    def f(x):
        return x*func(x/p)/p**2 * sp.jv(n, x) 
    
    for i in range(L):
       value += simpson(f, zeros[i], zeros[i+1], 1000)
    return value 

def f1(x): 
    return x**2
p = 1
n = 1 
L = 3000

p_value = Hankel_transform(f1, p, n, L)
print(p_value)
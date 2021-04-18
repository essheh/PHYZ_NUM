# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 21:36:40 2021

@author: franc
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator

def leibniz(n):
    value = 0
    part_sum = []
    for i in range(n):
        value += (-1) ** i * 4/(2*i+1)
        part_sum.append(value)
    
    
    return part_sum


def ln3(n):
    value = 0
    part_sum = []
    for i in range(n):
        value += ((-1)**(i)) * 2**(i+1) / (i+1)
        part_sum.append(value)
    
    
    return part_sum

part_sum = np.array(leibniz(13))
part_sum_ln = np.array(ln3(20))

def Shanks_transform(seq):
    """
    Return the Shanks transformation of a sequence of partial sums

    Parameters
    ----------
    seq : list 
        partial sum of the serie 

    Returns
    -------
    Shanks transformation values

    """
    part_shanks = []
    for i in range(1, len(seq)-1):
        value = seq[i+1] - (seq[i+1]- seq[i])**2/((seq[i+1]-seq[i])-(seq[i]-seq[i-1])) 
        part_shanks.append(value)
    return part_shanks

part_shanks_1 = np.array(Shanks_transform(part_sum))


part_shanks_2 = np.array(Shanks_transform(part_shanks_1))


part_shanks_3 = np.array(Shanks_transform(part_shanks_2))



def error_abs(approx, real) :
    return abs(approx - real)

plt.rcParams.update({'font.size': 14})
params = {'legend.fontsize': 9,
          'legend.handlelength': 2}
plt.rcParams.update(params)
plt.yscale("log")
plt.plot(np.arange(len(part_sum))+1, error_abs(part_sum, np.pi), label = "$A_n$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_1))+3, error_abs(part_shanks_1, np.pi), label = "$S(A_n)$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_2))+5, error_abs(part_shanks_2, np.pi), label = "$S^2(A_n)$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_3))+7, error_abs(part_shanks_3, np.pi), label = "$S^3(A_n)$", marker='.', linestyle = '-')
plt.ylabel('Erreur absolue [-]')
plt.xlabel('n')
plt.tight_layout()
plt.legend()
ax=plt.axes()
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
plt.savefig('figures/Shanks_lebniz', dpi=600)
plt.show()

part_shanks_ln_1 = np.array(Shanks_transform(part_sum_ln))
part_shanks_ln_2 = np.array(Shanks_transform(part_shanks_ln_1))
part_shanks_ln_3 = np.array(Shanks_transform(part_shanks_ln_2))

plt.yscale("log")
plt.plot(np.arange(len(part_sum_ln))+1, error_abs(part_sum_ln, np.log(3.0)), label = "$A_n$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_ln_1))+3, error_abs(part_shanks_ln_1, np.log(3.0)), label = "$S(A_n)$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_ln_2))+5, error_abs(part_shanks_ln_2, np.log(3.0)), label = "$S^2(A_n)$", marker='.', linestyle = '-')
plt.plot(np.arange(len(part_shanks_ln_3))+7, error_abs(part_shanks_ln_3, np.log(3.0)), label = "$S^3(A_n)$", marker='.', linestyle = '-')
plt.ylabel('Erreur absolue [-]')
plt.xlabel('n')
plt.tight_layout()
plt.legend()
ax=plt.axes()
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
plt.savefig('figures/Shanks_ln3', dpi=600)
plt.show()

def epsilon(terms):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    n = len(terms)
    e = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        e[i, 1] = terms[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    sumation = e[:, 1:n + 1:2]
    return sumation[-1,-1]

def leibniz(n):
    value = 0
    part_sum = []
    for i in range(n):
        value += (-1) ** i * 4/(2*i+1)
        part_sum.append(value)
    
    
    return part_sum


def ln3(n):
    value = 0
    part_sum = []
    for i in range(n):
        value += ((-1)**(i)) * 2**(i+1) / (i+1)
        part_sum.append(value)
    
    
    return part_sum
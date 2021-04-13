#Module
import numpy as np
import scipy.special as sp
import scipy.integrate as integ

def H(n,p,g):                                            # H_n(p;g)
        
    def f(x):                                            # f(x) = xg(x/p)/p^2
        return x * g(x/p) / p ** 2   
    
    func = f(x) * sp.jv(n,x)                             # f(x)*J_n(x)
    
    def zero(n, l):                                      # j_{n,l}
        return sp.jn_zeros(n, l)
    
    
    def integrale(func, n, l):                           # Intégrale
        return integ.quad(func,zero(n,l),zero(n,l+1))
    
    s = 0.0
    for l in range(0, 10000):                            # Somme des intégrale
        s += integrale(func,n,l)
    return s
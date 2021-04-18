import numpy as np
import mpmath
from mpmath import *
import matplotlib.pyplot as plt
import timeit
def leibniz(n):
    terms_leibniz = []
    part_sum_leibniz = []
    t_sum = 0
    for i in range(n):
        term = (-1) ** i * 4/(2*i+1)
        terms_leibniz.append(term)
        t_sum = t_sum + term
        part_sum_leibniz.append(t_sum)
  
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_pi_pct = 100*abs(np.pi - part_sum_leibniz[n-1]) / np.pi    
    
    return terms_leibniz, part_sum_leibniz, err_pi_pct


def ln3(n):
    terms_ln3 = []
    part_sum_ln3 = []
    t_sum = 0
    for i in range(n):
        term = ((-1)**(i)) * 2**(i+1) / (i+1)
        terms_ln3.append(term)
        t_sum = t_sum + term
        part_sum_ln3.append(t_sum)
        
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_ln3_pct = 100*abs(np.log(3) - part_sum_ln3[n-1]) / np.log(3)    
    
    return terms_ln3, part_sum_ln3, err_ln3_pct
    

def epsilon(terms):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    n = len(terms)
    e = np.zeros((n + 1, n + 1))
    for i in range(1, n + 1):
        e[i, 1] = terms[i - 1]#All the partial sum in the firts col

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    mat = e[:, 1:n + 1:2]
    if n>4:
        err_approx = abs(mat[-1,-3]-mat[-1,-1])
    else:
        err_approx = float(nan)
    return mat,err_approx,mat[-1,-1]

def compute(n0,nmax,func):
    """
    The function returns all the important output compute with the partial sumation and the epsilon algorithm..
    param:  n0: Lower bound nonzero n terms. 
            nmax:Upper bound of n terms.
            func: Summation function.
    return: n: number of terms
            terms: terms of the serie.
            sp, sp_err: Partial sum and partial sum absolute error.
            epsi, epsi_err: Sum compute with epsilon and the absolute error.
            epsi_err_approx: Approximative error 
    """
    n = []
    sp = []
    sp_err =[]
    terms = []
    epsi = []
    epsi_err =[]
    epsi_err_approx = []
    for N in range(n0,nmax):    
        term, s, err= func(N)
        terms.append(term) #List of N terms of the serie
        sp_err.append(err) #List of N absolute error of the partial sum
        sp.append(s[N-1]) #List of N partial sum of the serie
        if N >3:
            mat,err_approx,best= epsilon(s)
        else:
            err_approx = float(nan)
            best = float(nan)
        epsi_err_approx.append(err_approx) #List of N approx error
        epsi.append(best)#List of N best acurate value with epsilon algo
        a = 100*abs(np.pi - best) / np.pi 
        epsi_err.append(a)#List of N absolute error with epsilon algo
        n.append(N)
    return n,terms[-1],[sp,sp_err],[epsi,epsi_err],epsi_err_approx,mat
    

#Comparaison mpmat
K = 20
S = [4*sum(mpf(-1)**n/(2*n+1) for n in range(m)) for m in range(1,K)]
T = shanks(S)
test = compute(3, K, leibniz)
nprint(T[-1][-1])
print(test[3][0][-1])

t1 = timeit.timeit('[epsilon(T[1])]', number=1,globals=globals())
t2 = timeit.timeit('[shanks(S)]', number=1,globals=globals())
print(t1,t2)

# =============================================================================
# Fig
# =============================================================================
#FIGURES SOMMES PARTIELLES
Leb = compute(1,100,leibniz)
ln = compute(1,100,ln3)

plt.rcParams.update({'font.size': 17})
params = {'legend.fontsize': 10,
          'legend.handlelength': 2}
plt.rcParams.update(params)

plt.plot(Leb[0],Leb[2][0])
plt.ylabel('$A_n$')
plt.xlabel('n')
plt.ylim(2.5,4)
plt.xlim(0,100)
#plt.savefig('figures/leibniz_sum', dpi=600,bbox_inches='tight')
#plt.show()


plt.plot(ln[0],ln[2][0])
plt.yscale("symlog")
plt.ylabel('$A_n$')
plt.xlabel('n')
plt.yticks([-10e29,-10e19,-10e9,0,10e9,10e19,10e29])
#plt.savefig('figures/ln3_sum', dpi=600,bbox_inches='tight')
#plt.show()

#Erreur approximative
Leb2 = compute(1,50,leibniz)
ln2 = compute(1,50,ln3)
N = np.array(Leb2[0])*2
plt.plot(N, Leb2[4],label='Leibniz')
plt.plot(N, ln2[4],label='Ln3')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel("Erreur approximative")
plt.legend()
plt.xlim(0,50)
#plt.savefig('figures/err_approx', dpi=600,bbox_inches='tight')
#plt.show()

#Erreur relatives

plt.plot(N,ln2[3][1],'.', label = 'Algorithme $\epsilon$')
plt.plot(N,ln2[2][1],'.',label = 'Sommes partielles')
plt.legend()
plt.ylim(10e-6,10e8)
plt.xlim(0,50)
plt.ylabel('Erreur relative [%]')
plt.xlabel('N')
plt.yscale("log")
#plt.savefig('figures/ln3_err', dpi=600,bbox_inches='tight')
#plt.show()

plt.plot(N,Leb2[3][1],'.', label = 'Algorithme $\epsilon$')
plt.plot(N,Leb2[2][1],'.',label = 'Sommes partielles')
plt.ylabel('Erreur relative [%]')
plt.xlabel('N')
plt.yscale("log")
plt.xlim(0,50)
plt.tight_layout()
plt.legend()
#plt.savefig('figures/leibniz_err', dpi=600,bbox_inches='tight')
#plt.show()




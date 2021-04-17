import numpy as np
import mpmath
from mpmath import *
import matplotlib.pyplot as plt

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
    k = len(terms)
    n = int((np.array(k)-1)/2)
    e = np.zeros((n + 1, n + 1))
    for i in range(1, n + 1):
        e[i, 1] = terms[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])
    mat_epsi = e[:, 1:n + 1:2]
    if n>5:
        err_approx = abs(mat_epsi[-1,-3]-mat_epsi[-1,-1])
    else:
        err_approx = nan
    return mat_epsi,err_approx,mat_epsi[-1,-1]

k = 20
terms_leibniz, part_sum_leibniz, err_pi = leibniz(k)
mat,err,epsi = epsilon(part_sum_leibniz)



sum_leibniz = []
sum_ln3 = []
sum_leibniz_err = []
sum_ln3_err = []

epsi_leibniz = []
epsi_leibniz_err = []
err_aprox_leibniz = []
terms_leibniz = []

epsi_ln3 = []
epsi_ln3_err = []
err_aprox_ln3 = []
terms_ln3 =[]

for N in range(3,k):    
    terms, part_sum_leibniz, err_pi = leibniz(N)
    terms_leibniz.append(terms)
    sum_leibniz_err.append(err_pi)
    sum_leibniz.append(part_sum_leibniz[N-1])
    
    mat,err,epsi = epsilon(part_sum_leibniz)
    err_aprox_leibniz.append(err)
    err = 100*abs(np.pi - epsi) / np.pi
    epsi_leibniz_err.append(err)
    epsi_leibniz.append(epsi)
    
    terms, part_sum_ln3, err_ln3 = ln3(N)
    terms_ln3.append(terms)
    sum_ln3_err.append(err_ln3)
    sum_ln3.append(part_sum_ln3[N-1])
    
    mat,err,epsi  = epsilon(part_sum_ln3)
    err_aprox_ln3.append(err)
    err = 100*abs(np.log(3) - epsi) / np.log(3)
    epsi_ln3_err.append(err)
    epsi_ln3.append(epsi)
    
f = list(range(1, k))
f2 = 2*np.array(f)+1

sum_leibniz=[]
sum_ln3 = []
for N in range(1,k):    
    terms, part_sum_leibniz, err_pi = leibniz(N)
    sum_leibniz.append(part_sum_leibniz[N-1])
    terms, part_sum_ln3, err_ln3 = ln3(N)
    sum_ln3.append(part_sum_ln3[N-1])
    
# plt.rcParams.update({'font.size': 18})
# params = {'legend.fontsize': 10,
#           'legend.handlelength': 2}
# plt.rcParams.update(params)
# plt.plot(f, err_aprox_leibniz,label='Leibniz')
# plt.plot(f, err_aprox_ln3,label='Ln3')
# plt.yscale('log')
# plt.xlabel('N')
# plt.ylabel("Erreur approximative")
# plt.legend()
# plt.xlim(0,50)
# plt.savefig('figures/err_approx', dpi=600,bbox_inches='tight')
# plt.show()

plt.rcParams.update({'font.size': 12})
params = {'legend.fontsize': 10,
          'legend.handlelength': 2}
plt.rcParams.update(params)
#plt.plot(f2,epsi_leibniz, label = 'epsilon')
plt.plot(f,sum_leibniz,label = 'somme')
plt.ylabel('Somme partielle')
plt.xlabel('N')
plt.ylim(2.5,4)
plt.xlim(0,20)
plt.tight_layout()
#plt.legend()
#plt.savefig('figures/leibniz_sum', dpi=600,bbox_inches='tight')
plt.show()


#plt.plot(f2,epsi_ln3, label = 'epsilon')
plt.plot(f,sum_ln3,label = 'somme')
plt.yscale("symlog")
plt.ylabel('Somme partielle')
plt.xlabel('N')
plt.tight_layout()
plt.xlim(0,20)
plt.ylim(-10e4,10e4)
#plt.legend()
plt.savefig('figures/ln3_sum', dpi=600,bbox_inches='tight')
plt.show()


# plt.plot(f2,epsi_ln3_err,'.', label = 'Algorithme $\epsilon$')
# plt.plot(f,sum_ln3_err,'.',label = 'Sommes partielles')
# plt.legend()
# plt.xlim(0,50)
# plt.ylabel('Erreur relative [%]')
# plt.xlabel('N')
# plt.tight_layout()
# plt.yscale("log")
# #plt.savefig('figures/ln3_err', dpi=600,bbox_inches='tight')
# plt.show()

# plt.plot(f2,epsi_leibniz_err,'.', label = 'Algorithme $\epsilon$')
# plt.plot(f,sum_leibniz_err,'.',label = 'Sommes partielles')
# plt.ylabel('Erreur relative [%]')
# plt.xlabel('N')
# plt.yscale("log")
# plt.xlim(0,50)
# plt.tight_layout()
# plt.legend()
# #plt.savefig('figures/leibniz_err', dpi=600,bbox_inches='tight')
# plt.show()




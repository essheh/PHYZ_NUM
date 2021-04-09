import numpy as np

def leibniz(n):
    terms_leibniz = []
    t_sum = 0
    for i in range(n):
        term = (-1) ** i * 4/(2*i+1)
        terms_leibniz.append(term)
        t_sum = t_sum + term
        part_sum_leibniz = t_sum
        
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_pi_pct = 100*abs(np.pi - part_sum_leibniz) / np.pi    
    
    return terms_leibniz, part_sum_leibniz, err_pi_pct

n = 100

terms_leibniz, part_sum_leibniz, err_pi = leibniz(n)

print(terms_leibniz)

print(part_sum_leibniz)
print(err_pi)
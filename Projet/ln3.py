import numpy as np


def ln3(n):
    terms_ln3 = []
    t_sum = 0
    for i in range(n):
        term = ((-1)**(i)) * 2**(i+1) / (i+1)
        terms_ln3.append(term)
        t_sum = t_sum + term
        part_sum_ln3 = t_sum
        
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_ln3_pct = 100*abs(np.log(3) - part_sum_ln3) / np.log(3)    
    
    return terms_ln3, part_sum_ln3, err_ln3_pct

n = 1000

terms_ln3, part_sum_ln3, err_ln3 = ln3(n)

#print(terms_ln3)


print(part_sum_ln3)

print(err_ln3)
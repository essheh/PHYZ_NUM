import numpy as np

s=0
ps=[]
n=10
for k in range(0, n):
    term = 4 * np.power(-1, k) / (2 * k + 1)
    s += term
    ps.append(s)
    
def epsilon(S):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    n = len(S)
    e = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        e[i, 1] = S[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    sum = e[:, 1:n + 1:2]
    return sum

sum = epsilon(ps)
print(sum[-1, -1])



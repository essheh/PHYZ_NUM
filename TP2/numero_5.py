import numpy as np
from numpy import exp, linspace
from pylab import plot, show, xlabel, ylabel, grid
np.seterr(divide = 'ignore')

#Constantes

R0_max = 5.7
accuracy = 1e-6
steps = 0.01

def p_graph(R0_max, accuracy, steps):
    

#--Calcul itératif du pourcentage de la population qui aura été infecté au terme de la pandémie

    #--Paramètres:
    #   R0_max: Nombre de reproduction de base maximale en supposant que le fluctuation commence à 0
    #   accuracy: précision de l'erreur de p
    #   steps: Bonds par itérations

    #--Retourne un graphique de p en fonction de R0

    points = int(R0_max / steps)


# Liste pour le graphique

    y = []
    temp = linspace(0, R0_max, points)

# R0 loop

    for R0 in temp:
        p1 = 1.0
        error = 1.0
    
    # Loop jusqu'à l'erreur établie
    
        while error > accuracy:
            p1,p2 = 1 - exp(-R0*p1),p1
            error = abs((p1-p2)/(1 - 1 / (R0 * exp(-R0 * p1))))
        
        # Ajout des valeurs convergées et transformation en pourcentage
        
        y.append(100 * p1)
    
#Graphique    
    
    plot(temp,y)
    xlabel("Nombre de reproduction de base ($R_0$)")
    ylabel("Pourcentage de la population qui aura \n été infecté au terme de la pandémie (p) [%]")
    grid()
    show()
        
p_graph(5.7,1e-6,0.01)
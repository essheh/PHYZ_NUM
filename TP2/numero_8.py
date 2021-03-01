def over_relax(R0, omega, accuracy):
    #--Calcul du nombre d'itération du pourcentage de la population qui aura été infecté au terme de la pandémie
    #--Selon la méthode de relaxation accélérée

    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   omega: paramètre de relaxation accélérée
    #   accuracy: précision de l'erreur de p

    #--Retourne le nombre d'itération selon la méthode de relaxation accélérée
    
    iterations = 1

    def f(p):
        return 1 - exp(- R0 * p)

    def f_derive(p):
        return R0 * exp(- R0 * p)

    def error(p1, p2):
        return (p1 - p2) / (1 - 1 / ((1 + omega) * f_derive(p1) - omega))

    p1 = 0.5  # valeur de départ
    p2 = (1 + omega) * f(p1) - omega * p1
    while abs(error(p1, p2)) > accuracy:
        p1, p2 = p2, (1 + omega) * f(p2) - omega * p2
        iterations += 1
    return iterations

def over_graph(R0,accuracy):
    #-- Calcul le graphique du nombre d'itération de la méthode de relaxation accélérée selon le paramètre omega partant 
    #-- de 0.5 en allant jusqu'à 2 en bond de 0.1
    
    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   accuracy: précision de l'erreur de p

    #--Retourne un graphique du nombre d'itération selon le paramètre de relaxation accélérée
    
    
    y = []
    temp = linspace(0.5,2,16)
    
    for omega in temp:
        iterations = over_relax(R0,omega, accuracy)
        
        y.append(iterations)
        
    #Graphique    
    
    plot(temp,y)
    xlabel("Valeur de $\omega$")
    ylabel("Nombre d'itérations")
    grid()
    show()
    
over_graph(2,1e-6)

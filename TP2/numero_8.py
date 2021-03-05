def over_relax(R0, omega,start, accuracy):
    #--Calcul du nombre d'itération du pourcentage de la population qui aura été infecté au terme de la pandémie
    #--Selon la méthode de relaxation accélérée

    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   omega: paramètre de relaxation accélérée
    #   start: valeur de départ pour l'itération allant de 0.1 à 1.4
    #   accuracy: précision de l'erreur de p

    #--Retourne le nombre d'itération, la valeur finale de p selon la méthode de relaxation accélérée,
    #--la valeur de R0 ainsi que la valeur de départ
    
    iteration = 1

    def f(p):
        return 1 - exp(- R0 * p)

    def f_derive(p):
        return R0 * exp(- R0 * p)

    def error(p1, p2):
        return (p1 - p2) / (1 - 1 / ((1 + omega) * f_derive(p1) - omega))

    p1 = start  # valeur de départ
    p2 = (1 + omega) * f(p1) - omega * p1
    while abs(error(p1, p2)) > accuracy:
        p1, p2 = p2, (1 + omega) * f(p2) - omega * p2
        iteration += 1
    return iteration,p2,R0,start

def over_graph(R0,start,accuracy):
    #-- Calcul le graphique du nombre d'itération de la méthode de relaxation accélérée selon le paramètre omega partant 
    #-- de 0.5 en allant jusqu'à 2 en bond de 0.1
    
    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   start: valeur de départ pour l'itération allant de 0.1 à 1.4
    #   accuracy: précision de l'erreur de p

    #--Retourne un graphique du nombre d'itération selon le paramètre de relaxation accélérée
    
    
    y = []
    temp = linspace(0.5,2,16)
    
    for omega in temp:
        iterations = over_relax(R0,omega,start, accuracy)[0]
        
        y.append(iterations)
        
    #Graphique    
    
    plot(temp,y)
    xlabel("Valeur de $\omega$")
    ylabel("Nombre d'itérations")
    grid()
    show()
    
over_graph(2,1,1e-6)
print("Selon la courbe ci-dessus, pour un nombre de reproduction de base (R_0) de {:} et une valeur de départ de {:}, la valeur \noptimale de omega est 0.5 ou 0.7, afin d'obtenir une valeur de p de {:.2f}% en seulement {:} itérations".format(over_relax(2, 0.5,1, 1e-6)[2],over_relax(2, 0.5,1, 1e-6)[3],over_relax(2, 0.5,1, 1e-6)[1]*100,over_relax(2, 0.5,1, 1e-6)[0]))

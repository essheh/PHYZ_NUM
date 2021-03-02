def p_nb_iter(R0,start,accuracy):
    
    #--Calcul du nombre d'itération du pourcentage de la population qui aura été infecté au terme de la pandémie
    #--Selon la méthode de relaxation

    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   start: Valeur de départ pour l'itérations allant de 0.1 à 1.4
    #   accuracy: précision de l'erreur de p

    #--Retourne le nombre d'itération et la solution convergente à la précision préétablie
    
    iteration = 1
    def f(p):
        return 1 - exp(- R0 * p)

    def error(p1, p2):
        return (p1 - p2) / (1 - 1 / (R0 * exp(-R0 * p1)))

    p1 = start # valeur de départ
    p2 = f(p1)
    while(abs(error(p1, p2)) > accuracy):
        p1, p2 = p2, f(p2)
        iteration += 1
    print("Nombre d'itérations : ", iteration)
    return p2

print(p_nb_iter(2,1,1e-6))
def p_nb_iter(R0,start,accuracy):
    
    #--Calcul du nombre d'itération du pourcentage de la population qui aura été infecté au terme de la pandémie
    #--Selon la méthode de relaxation

    #--Paramètres:
    #   R0: Nombre de reproduction de base 
    #   start: Valeur de départ pour l'itérations allant de 0.1 à 1.4
    #   accuracy: précision de l'erreur de p

    #--Retourne le nombre d'itération, la solution convergente à la précision préétablie selon la méthode de relaxation,
    #--la valeur de R0 et la valeur de départ pour en faire des variables globales
    
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
    return iteration,p2,R0,start


print("Pour un nombre de reproduction de base (R_0) de {:} et une valeur de départ de {:}, on obtient un pourcentage \nde population (p) de {:.2f} %, en {:} itérations.".format(p_nb_iter(2,0.5,1e-6)[2],p_nb_iter(2,0.5,1e-6)[3],p_nb_iter(2,0.5,1e-6)[1]*100,p_nb_iter(2,0.5,1e-6)[0]))
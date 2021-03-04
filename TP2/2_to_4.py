import matplotlib.pyplot as plt
from scipy import integrate
import numpy as np
from math import sqrt, pi 

def LJP(x):
    """
    The function return the adimensional Leonard-Jones Potential values
    
    param 1 x: x is the adimensional distance value 
    return: adimensional Leonard-Jones Potential value 
    """ 
    return 4*((1/x)**12 - (1/x)**6)
    
x = np.linspace(0.9,4,1000)

plt.plot(x, LJP(x))
#plt.show()

def turning_points(epsilon):
    """
    The function gives the turning points of the adimensional Leonard-Jones Potential
    for different energy values
    
    param 1 epsilon: epsilon is the adimensional energy value of which the turning points are searched 
    return: The two turning point of the function for the energy value  
    """ 
   
    first_tp = 2**(1/6) * (sqrt(epsilon+1)/epsilon-1/epsilon)**(1/6)
    second_tp = 2**(1/6) * (-((sqrt(epsilon+1)+1)/epsilon))**(1/6)
    return first_tp, second_tp
    
def LJP_discrete_energy(n, gamma):
    """
    The function gives the n first adimensional discrete energy values using the WKB approximation the solve the Schrödinger equation
    for the Lennard-Jones potential, using the bissection method 
    
    param 1 n: n is the number of energy levels to find, it must be a positive integer   
    param 2 gamma: gamma is a constant that varies between each molecule, it must have a positive value 
    return: The n first adimensional discrete energy values   
    """ 
    energies = [] #The list of the n first discrete energy values 
    accuracy = 10E-08 #The accuracy 
    def f(epsilon, n): # function that return the value of the function in the equation (3)  
        func = lambda x: (epsilon - LJP(x))**0.5 # The function to integrate, i.e., (E - V(x))**(1/2) 
        integral = lambda epsilon: integrate.quad(func, a = turning_points(epsilon)[0], b = turning_points(epsilon)[1])[0] #the integral at the turning_points
        return integral(epsilon) - (n+0.5)*pi/gamma
    
    def midpoint(x1, x2): 
        return (x1 + x2) / 2

    def have_same_sign(x1, x2): #function used for the bissection method 
        if x1 < 0 and x2 < 0 or x1 > 0 and x2 > 0:
            return True
        else:
            return False
    
    for i in range(n): # we repeat the loop for the n energies
        
        #Initial values of the bissection method 
        x1 = -0.9999999
        x2 = -0.0000001
        
        # From here we start the bissection method 
        
        while abs(x1 - x2) > accuracy: 
            x = midpoint(x1, x2)
            if have_same_sign(f(x1, i), f(x, i)):
                x1 = x
            elif have_same_sign(f(x, i), f(x2, i)):
                x2 = x
            elif abs(x) < accuracy:
                return x
        
        energies.append(midpoint(x1, x2))
    return energies 

#print(LJP_discrete_energy(20, 150))

def plot_energy_levels(n, gamma): 
    x_LJP = np.linspace(1.001,1.9,1000)
    for i in range(n):
        epsilon = LJP_discrete_energy(n, gamma)[i]
        #Plot the energy levels 
        plt.plot(np.linspace(turning_points(epsilon)[0], turning_points(epsilon)[1], 1000), np.full(1000, epsilon) , color = "black", linestyle = "dashed") 
        plt.plot(turning_points(epsilon)[0], epsilon, color = "red",  marker='o') #Plot the first TP
        plt.plot(turning_points(epsilon)[1], epsilon, color = "red",  marker='o') #Plot the second TP 
    plt.plot(x_LJP, LJP(x_LJP), color = "black") #Plot the Lennard-Jones potential
    plt.show()

#plot_energy_levels(20, 150)

def LJP_discrete_energy_Secant(n, gamma):
    """
    The function gives the n first adimensional discrete energy values using the WKB approximation the solve the Schrödinger equation
    for the Lennard-Jones potential, using the secant method 
    
    param 1 n: n is the number of energy levels to find, it must be a positive integer   
    param 2 gamma: gamma is a constant that varies between each molecule, it must have a positive value 
    return: The n first adimensional discrete energy values   
    """ 
    energies = [] #The list of the n first discrete energy values 
    def f(epsilon, n): # function that return the value of the function in the equation (3)  
        func = lambda x: (epsilon - LJP(x))**0.5 # The function to integrate, i.e., (E - V(x))**(1/2) 
        integral = lambda epsilon: integrate.quad(func, a = turning_points(epsilon)[0], b = turning_points(epsilon)[1])[0] #the integral at the turning_points
        return integral(epsilon) - (n+0.5)*pi/gamma
    
    for i in range(n): # we repeat the loop for the n energies
        
        #Initial values of the secant method 
        x1 = -0.9999999
        x2 = -0.0000001
        
        # From here we start the secant method 
        for j in range(8):
            x3 = x2 - f(x2, i) * (x2 - x1) / float(f(x2,i) - f(x1, i))
            x1, x2 = x2, x3
        energies.append(x3)
    
    return energies 
print(LJP_discrete_energy(20, 150))
print(LJP_discrete_energy_Secant(20, 150))
        
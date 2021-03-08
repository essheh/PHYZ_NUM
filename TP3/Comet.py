from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt

#Constants  

#M = const.M_sun #Mass of the sun [kg]
#G = const.G #Gravitationnal constant [m^3/(kg s^2)]
#G_year = G*(365.2422*24*60*60)**2 #Gravitationnal constant [m^3/(kg year^2)] 

G_year = 66465.32132945514
M = 1.988409870698051e+30

def f(r, t): #f returns the vector of the derivatives of each coordinates 
    x = r[0]
    y = r[1]
    v_x = r[2]
    v_y = r[3]
    dist = (x**2 + y**2)**(0.5)
    return np.array([v_x, v_y, -G_year*M*x/dist**3, -G_year*M*y/dist**3], float)

def comet_orbit(x, y, v_x, v_y, t, h):
    """
    The function return the x and y position at time t of the orbit of a comet using the 4th order Runge-Kutta using Newton's second law.    
    
    param 1 x: x is the distance in the x plane between the sun and the comet à time t = 0 [m] 
    param 2 y: y is the distance in the y plane between the sun and the comet à time t = 0 [m] 
    param 3 v_x: v_x is the speed of the comet in the x direction at time t = 0 [m/year] 
    param 4 v_y: v_y is the speed of the comet in the y direction at time t = 0 [m/year] 
    param 5 t: t is the the time duration of the simulation [year]
    param 6 h: h is the time step between each time value 
    return: x_points, y_points and time two array of the x and y position for each time t and one of the time duration with step h  
    """ 
    t_points = np.arange(0, t, h) #The array of time values 
    x_points = [] #The list of x positions values
    y_points = [] #The list of y positions values
    r = np.array([x, y, v_x, v_y], float) #The vector of all coordinates  
  
    #The Runge-Kutta method starts from here 
    for t in t_points:
        x_points.append(r[0])
        y_points.append(r[1])
        
        k_1 = h*f(r, t)
        k_2 = h*f(r + 0.5*k_1, t + 0.5*h)
        k_3 = h*f(r + 0.5*k_2, t + 0.5*h)
        k_4 = h*f(r + k_3, t + h)
        r += (k_1 + 2*k_2 + 2*k_3 + k_4)/6 
    
    return x_points, y_points
 
x = 4*10**12 #[m]
y = 0 
v_x = 0
v_y = 500*(365.2422*24*60*60) # [m/year]
t = 80 #[year]
h = 0.0008

x_points, y_points = comet_orbit(x, y, v_x, v_y, t, h)
plt.plot(x_points, y_points)
plt.show()






from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt

#Constants  

#M = const.M_sun #Mass of the sun [kg]
#G = const.G #Gravitationnal constant [m^3/(kg s^2)]
#G_year = G*(365.2422*24*60*60)**2 #Gravitationnal constant [m^3/(kg year^2)] 

G_year = 66465.32132945514
M = 1.988409870698051e+30

def f(r): #f returns the vector of the derivatives of each coordinates 
    x = r[0]
    y = r[1]
    v_x = r[2]
    v_y = r[3]
    dist = (x**2 + y**2)**(0.5)
    return np.array([v_x, v_y, -G_year*M*x/dist**3, -G_year*M*y/dist**3], float) 

def Runge_Kutta(f, r, h):
    """ 
    The function return the vector r of the updated coordinates using the Runge-Kutta method

    param 1 f: f is a vector of function
    param 2 r: r is the vector of the current coordinates
    param 3 h: h is the step between each coordinates 
    return: r the updated coordinates vector 
    """
    k_1 = h*f(r)
    k_2 = h*f(r + 0.5*k_1)
    k_3 = h*f(r + 0.5*k_2)
    k_4 = h*f(r + k_3)
    r += (k_1 + 2*k_2 + 2*k_3 + k_4)/6 
    return r
    

def orbit_constant_step(x, y, v_x, v_y, t, h):
    """
    The function return the x and y position at time t of the orbit of a comet using the 4th order Runge-Kutta to solve Newton's law of universal
    gravitation with constant time step h.    
    
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
    r = np.array([x, y, v_x, v_y]) #The vector of all coordinates  
  
    #The Runge-Kutta method starts from here 
    for t in t_points:
        x_points.append(r[0])
        y_points.append(r[1])   
        r = Runge_Kutta(f, r, h) 
    
    return x_points, y_points
 
x = 4*10**12 #[m]
y = 0 
v_x = 0
v_y = 1000*(365.2422*24*60*60) # [m/year]
t = 40 #[year]
h = 0.0008

x_points, y_points = orbit_constant_step(x, y, v_x, v_y, t, h)
plt.plot(x_points, y_points)
plt.plot(0, 0, color = 'yellow', marker='o')
plt.show()


def time_step(x, y, v_x, v_y, t, h, delta):
    """
    The function return the x and y position at time t of the orbit of a comet using the 4th order Runge-Kutta to solve Newton's law of universal
    gravitation with variable time step.    
    
    param 1 x: x is the distance in the x plane between the sun and the comet à time t = 0 [m] 
    param 2 y: y is the distance in the y plane between the sun and the comet à time t = 0 [m] 
    param 3 v_x: v_x is the speed of the comet in the x direction at time t = 0 [m/year] 
    param 4 v_y: v_y is the speed of the comet in the y direction at time t = 0 [m/year] 
    param 5 t: t is the the time duration of the simulation [year]
    param 6 h: h is the time step between each time value 
    param 7 delta: delta is the maximal accepted error  
    return: x_points, y_points and time two array of the x and y position for each time t and one of the time duration with step h  
    """
    
    t_points = np.arange(0, t, h) #The array of time values 
    x_points = [] #The list of x positions values
    y_points = [] #The list of y positions values
    r = np.array([x, y, v_x, v_y]) #The vector of all coordinates 
    
    step_1 = Runge_Kutta(f, r, h)
    step_2 = Runge_Kutta(f, r + step_1, h)
    total_step = step_1 + step_2
    
    double_step = Runge_Kutta(f, r, 2*h)
    
    x_error = total_step[0] - double_step[0]
    y_error = total_step[1] - double_step[1]
    error = np.sqrt(x_error**2 + y_error**2) 
    rho = 30*h*delta/(error)
    
    if rho > 1: 
            new_h = h*rho**(1/4)
    else: 
        
    
delta = 1000 #[m]    




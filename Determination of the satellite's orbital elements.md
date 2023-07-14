
1. Importing the required stuff
~~~python
import sympy as sp  
from sympy import symbols, solve  
from sympy.solvers import nonlinsolve  
import numpy as np  
from datetime import datetime, timedelta  
import matplotlib.pyplot as plt  
from matplotlib import animation  
from scipy.integrate import odeint
~~~

2. Solving the first part of this project
~~~python
def Task1 (x,y,z,Vx,Vy,Vz,M):

	# Calculates the magnitudes of position and velocity vectors
    mr = sp.sqrt(x ** 2 + y ** 2 + z ** 2)  
    mv = sp.sqrt(Vx ** 2 + Vy ** 2 + Vz ** 2)  
    
    # Converts the magnitudes to numerical values
    r = sp.N(mr)  
    v = sp.N(mv)  

	# Calculates the dot product between the position and velocity vectors
    sg = x * Vx + y * Vy + z * Vz  

	# Calculates the angle between the vectors using the arccosine function
    gm = sp.acos(sg / (r * v))  

	# Solves a system of equations using the nonlinsolve function, resulting in             values for `rp` and `vp`
    rp, vp = symbols('rp, vp', real=True)  
    p = list(nonlinsolve([r * v * sp.sin(gm) - rp * vp, vp ** 2 - v ** 2 - 2 * M      *     (1 / rp - 1 / r)], [rp, vp]))  

	# Assigns the values of `rp` and `vp` to variables `A` and `B`
    A = p[0]  
    B = p[1]  

	# Checks conditions and assigns new values to `p` based on the conditions
    if sp.N((r * v ** 2 / M - 1),1,chop=True)==0:  
        p[0]=r  
        p[1]=v  
    elif A[0] > B[0] or A[0]<0:  
        p = B  
    else:  
        p = A  

	# Calculates the eccentricity `e` using `p` values
    e = p[0] * p[1] ** 2 / M - 1  
    e = sp.N(e, 5)  

	# Calculates the true anomaly `nu`
    nu = sp.atan((r * v ** 2 / M) * sp.sin(gm) * sp.cos(gm) / (r * v ** 2 / M *             sp.sin(gm) ** 2 - 1))  

	# Calculates the inclination `i`
    sg = v * r * sp.sin(gm)  
    Cz = x * Vy - Vx * y  
    Cy = z * Vx - x * Vz  
    Cx = y * Vz - z * Vy  
    ir = sp.acos(Cz / sg)  
    i = ir * 180 / sp.pi  
    i = sp.N(i, 8)  

	# Calculates various forces components `Fx`, `Fy`, and `Fz`
    Fx=Vy*Cz-Vz*Cy-M*x/r  
    Fy=Vz*Cx-Vx*Cz-M*y/r  
    Fz=Vx*Cy-Vy*Cx-M*z/r  

	# Calculates the resultant force `F`
    F=sp.sqrt(Fx**2+Fy**2+Fz**2)  

	# Calculates unit vectors `Enx` and `Eny`
    Enx=-Cy/(sp.sqrt(Cx**2+Cy**2))  
    Eny=Cx/(sp.sqrt(Cx**2+Cy**2))  

	# Calculates the argument of the ascending node `omg`
    omgr=sp.acos(((Enx*Fx)+(Eny*Fy))/F)  
    omg = omgr * 180 / sp.pi  
    omg = sp.N(omg,8)  

	# Calculates the longitude of the ascending node `x`
    xr = sp.acos(-(Cy / sp.sqrt(Cx ** 2 + Cy ** 2)))  
    x = xr * 180 / sp.pi  
    x = sp.N(x,8)  

	# Calculates the angles `E0`, `EP`, `M0`, and `MP` related to the eccentric             anomaly
    E0 = sp.acos((e + sp.cos(nu)) / (1 + e * sp.cos(nu)))  
    EP = sp.acos((e + 1) / (1 + e))  
    M0 = E0 - e * sp.sin(E0)  
    MP = EP - e * sp.sin(EP)  

	# Checks conditions and assigns values to `Fp`, `Tp`, `a`, and `n`
    if sp.N(e,1)==1:  
        Fp='Nan'  
        Tp='Nan'  
        a=0  
        n=0  
    else:  
        a = p[0] / (1 - e)  
        #print("a", a)  
        #print("e", e)        Fp = a * (1 - e ** 2)  
        Fp = sp.N(Fp, 10)  
        n = sp.sqrt(M / a ** 3)  
        Tp = (MP - M0) / n  
        if Tp<=0:  
            Tp = sp.N(Tp,8)  
        else:  
            Tp = (M0 - MP) / n  
            Tp = sp.N(Tp, 8)  
        
    print('Focal Parameter =', Fp, "\nEccentricity =", e, "\nInclination =", i,             '\nArgument of Periapsis =', omg,  
          "\nLongitude of Ascending Node =", x)  
    print('Time of passage of the periapsis t=', Tp, 'c', )  
    return (Fp,e,ir,omgr,xr,M0,n,a)
~~~


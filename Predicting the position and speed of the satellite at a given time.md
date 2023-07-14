
3. Solving the second part
~~~python
def task2 (Fp,e,i,omg,x,N,M0,n,a,M):  

	# Calculates the time `t1` in seconds based on the provided `N` (number of              hours)
    t1=N*60*60

	# Calculates the mean anomaly `M1` at time `t1` (product of mean motion (`n`))
    M1=M0+n*t1  

	# Iteratively solves Kepler's equation using the Newton-Raphson method to               find the eccentric anomaly (`E`)
    E=0  
    Epr=-1  

	# the loop continues until the absolute difference between `E` and its                  previous value (`Epr`) is less than 0.001
    while abs(E-Epr)>0.001:  
        Epr=E  
        E=e*sp.sin(E)+M1  

	# once the eccentric anomaly `E` is determined, the true anomaly (`Nu`) is              calculated using trigonometric formulas
    Nu=2*sp.atan(sp.sqrt((1+e)/(1-e))*sp.tan(E/2))  

	# Calculates the satellite's distance from the Earth's center
    r=(a*(1-e**2)/(1+e*sp.cos(Nu)))  

	# Computes the velocity components in the radial and tangential directions
    vt=(sp.sqrt(M/Fp)*(1+e*sp.cos(Nu)))  
    vr=(sp.sqrt(M/Fp)*e*sp.sin(Nu))  
    R=np.array([[r],[0],[0]])  
    V=np.array([[vr],[vt],[0]]) 

	# Constructs a rotation matrix `A` based on the orbital elements, which transforms      the position and velocity vectors from the orbital frame to the geocentric              equatorial coordinate system
    A=np.array([[(sp.cos(omg+Nu)*sp.cos(x)-sp.sin(omg+Nu)*sp.cos(i)*sp.sin(x)),(-sp.sin(omg+Nu)*sp.cos(x)-sp.cos(omg+Nu)*sp.cos(i)*sp.sin(x)),(sp.sin(i)*sp.sin(x))],\  
                [(sp.cos(omg+Nu)*sp.sin(x))+sp.sin(omg+Nu)*sp.cos(i)*sp.cos(x),(-sp.sin(omg+Nu)*sp.sin(x)+sp.cos(omg+Nu)*sp.cos(i)*sp.cos(x)),(-sp.sin(i)*sp.cos(x))],\  
                [(sp.sin(omg+Nu)*sp.sin(i)),(sp.cos(omg+Nu)*sp.sin(i)),(sp.cos(i))]])  

	# The position vector `R` and velocity vector `V` are transformed using the             rotation matrix `A` via matrix multiplication
    R=np.matmul(A,R)  
    V=np.matmul(A,V)  

	# The position and velocity vectors are then flattened, converted to lists, and         assigned to the variables `R` and `V`
    R=R.flatten().tolist()  
    V=V.flatten().tolist()  

	# The position coordinates (`x`, `y`, `z`) and velocity components (`Vx`, `Vy`,         `Vz`) are obtained by extracting the respective values from the `R` and `V` lists
    x,y,z=sp.N(R[0]),sp.N(R[1]),sp.N(R[2])  
    Vx,Vy,Vz=sp.N(V[0]),sp.N(V[1]),sp.N(V[2])  
    print('After ', N, ' hour(s)')  
    print('Coordinates: х- ',x,' y- ',y,' z- ',z)  
    print('Velocities: Vx- ',Vx,' Vy- ',Vy,' Vz- ',Vz)  

	# Returns the eccentric anomaly
    return E
~~~

4. Function to calculate the topographic coordinates (longitude and latitude) of the satellite based on its geocentric equatorial coordinates and other parameters
~~~python
def Topographic_coordinates(R, om_e, dt, lam0, R0):  

	# Calculates the magnitude of the position vector
    r = sp.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])

	# Computes the difference in longitudes between the satellite's position and the        reference position (`R` and `R0`)
    dlam = sp.atan2(R[0], R[1]) - sp.atan2(R0[0], R0[1])  

	# The initial longitude `lam0` is adjusted by subtracting the difference in             longitudes (`dlam`) and the Earth's rotation rate multiplied by the time difference     (`om_e * dt`)
    lam = lam0 + dlam - om_e * dt  

	# The calculated longitude `lam` is then adjusted to ensure it falls within the         range of -π to π
    while lam > sp.pi: lam -= 2 * sp.pi  
    while lam < -sp.pi: lam += 2 * sp.pi  

	# Calculates the lattitude
    fi = sp.asin(R[2] / r)  
    return lam, fi
~~~

5. Constructing a 2D plot of the satellite's trajectory
~~~python
def task3 (Rr,lam0,a,M,e,i,omg,x,M0,n):  
    print(a)  

	# Calculates the orbital period of the satellites
    T = 2 * sp.pi * sp.sqrt(a ** 3 / M)  
    print("T", T)  

	# The variable `om_e` is set to the Earth's rotation rate in radians per second
    om_e=(360*sp.pi/180)/(23*60*60+56*60+4)  

	# The variables `t1` and `dt` are defined to set the initial time and time step, respectively, for iterating over the satellite's trajectory
    t1=3*60  
    dt=3*60  

	# Calculates the magnitude of the initial position vector
    r = sp.sqrt(Rr[0] * Rr[0] + Rr[1] * Rr[1] + Rr[2] * Rr[2])  

	# The latitude `fi` is determined by taking the arcsine of the ratio of the z-coordinate of the initial position vector to its magnitude
    fi=sp.asin(Rr[2] / r)  

	# Two empty lists, `X` and `Y`, are initialized to store the longitude and latitude coordinates of the satellite's trajectory
    X,Y=[sp.N(lam0*180/sp.pi,4)],[sp.N(fi*180/sp.pi,4)]  

	# Calculates the satellite's position at different points in time
    while t1<T*2:  
	    # Calculates the mean anomaly
        M1 = M0 + n * t1  
        E = 0  
        Epr = -1  
        
        # The eccentric anomaly `E` is calculated iteratively using an initial estimate of 0 and the equation `E = e * sin(E) + M1` until the difference between consecutive estimates is smaller than 0.001
        while abs(E - Epr) > 0.001:  
            Epr = E  
            E = e * sp.sin(E) + M1  
        Nu = 2 * sp.atan(sp.sqrt((1 + e) / (1 - e)) * sp.tan(E / 2))  
        r = (a * (1 - e ** 2) / (1 + e * sp.cos(Nu)))  
        R = np.array([[r], [0], [0]])  
        
        # The rotation matrix `A` is constructed using the values of `omg`, `x`, and `i` to transform the position vector from the orbital coordinate system to the Earth-fixed coordinate system
        A = np.array([[(sp.cos(omg + Nu) * sp.cos(x) - sp.sin(omg + Nu) * sp.cos(i) * sp.sin(x)),  
                       (-sp.sin(omg + Nu) * sp.cos(x) - sp.cos(omg + Nu) * sp.cos(i) * sp.sin(x)),  
                       (sp.sin(i) * sp.sin(x))], \  
                      [(sp.cos(omg + Nu) * sp.sin(x)) + sp.sin(omg + Nu) * sp.cos(i) * sp.cos(x),  
                       (-sp.sin(omg + Nu) * sp.sin(x) + sp.cos(omg + Nu) * sp.cos(i) * sp.cos(x)),  
                       (-sp.sin(i) * sp.cos(x))], \  
                      [(sp.sin(omg + Nu) * sp.sin(i)), (sp.cos(omg + Nu) * sp.sin(i)), (sp.cos(i))]])  
                      
        # The position vector `R` is obtained by multiplying the rotation matrix `A` with the vector `R` and flattening it to a 1D list.
        R = np.matmul(A, R)  
        R = R.flatten().tolist()  
        R= [sp.N(R[0]), sp.N(R[1]), sp.N(R[2])]  
        R0=Rr  
        
        # The longitude and latitude coordinates of the current position `R` are computed using the `Topographic_coordinates` function, with the longitude `lam0` and reference position vector `R0` as inputs
        Xx,Yy=Topographic_coordinates(R, om_e, dt, lam0, R0)  
        lam0 = Xx  
        Xx=sp.N(Xx*180/sp.pi,4)  
        Yy=sp.N(Yy*180/sp.pi,4)  
        X.append(Xx)  
        Y.append(Yy)  
        Rr = R  
        
		# The time `t1` is incremented by the time step `dt` to move to the next iteration
        t1=t1+dt  

	# After the loop, the lists `X` and `Y` are divided into three sets of coordinates based on longitude values: `X1`, `Y1`, `X2`, `Y2`, `X3`, `Y3`. This is done to handle cases where the trajectory crosses the boundary of -180° longitude
    X31=X  
    Y31=Y  
    X1,X2,Y1,Y2,X3,Y3=[],[],[],[],[],[]  
    a=0  
    for i in range (len(X)):  
  
        if -180<X[i]<-170:  
            a=1  
        if a==0:  
            X1.append(X[i])  
            Y1.append(Y[i])  
        else:  
            X2.append(X[i])  
            Y2.append(Y[i])  
    a=0  
    X=X2  
    Y=Y2  
    X2,Y2=[],[]  
    for i in range (len(X)):  
  
        if -180<X[i]<-170 and i>1:  
            a=1  
        if a==0:  
            X2.append(X[i])  
            Y2.append(Y[i])  
        else:  
            X3.append(X[i])  
            Y3.append(Y[i])  
    x = np.array(X1)  
    y = np.array(Y1)  
    x1 = np.array(X2)  
    y1 = np.array(Y2)  
    x2 = np.array(X3)  
    y2 = np.array(Y3)  
    img = plt.imread("Earth.jpg")  
    plt.imshow(img, extent=[-180, 180, -90, 90])  
    plt.grid()  
    plt.plot(x,y,'r-')  
    plt.plot(x1, y1, 'r-')  
    plt.plot(x2, y2, 'r-')  
    plt.show()  
    return X31,Y31
~~~

6. Constructing a 3D plot of the satellite's trajectory
~~~python
def task3_3d (ZVx,ZVy,ZVz,ZX,ZY,ZZ,M,):  

	# The function below represents the two-body problem. This function takes the state of the satellite and the time as inputs and calculates the derivatives of the state variables (position and velocity) based on the gravitational force from the Earth
    def model_2BP(state, t):  
        mu = M  
        x = state[0]  
        y = state[1]  
        z = state[2]  
        x_dot = state[3]  
        y_dot = state[4]  
        z_dot = state[5]  
        x_ddot = -mu * x / (x ** 2 + y ** 2 + z ** 2) ** (3 / 2)  
        y_ddot = -mu * y / (x ** 2 + y ** 2 + z ** 2) ** (3 / 2)  
        z_ddot = -mu * z / (x ** 2 + y ** 2 + z ** 2) ** (3 / 2)  
        dstate_dt = [x_dot, y_dot, z_dot, x_ddot, y_ddot, z_ddot]  
        return dstate_dt  

	# The initial position and velocity components (`X_0`, `Y_0`, `Z_0`, `VX_0`, `VY_0`, `VZ_0`) are assigned using the input parameters
    X_0 = ZX  
    Y_0 = ZY  
    Z_0 = ZZ  
    VX_0 = ZVx  
    VY_0 = ZVy  
    VZ_0 = ZVz  

	# The `state_0` variable is created as a list containing the initial state of the satellite
    state_0 = [X_0, Y_0, Z_0, VX_0, VY_0, VZ_0]  

	# The time `t` is defined as an array using `np.linspace()` to simulate the trajectory for a time period of 6 hours.
    t = np.linspace(0, 6 * 3600, 200)

	# The satellite's position components (`X_Sat`, `Y_Sat`, `Z_Sat`) are extracted from the `sol` array
    X_Sat = sol[:, 0]  
    Y_Sat = sol[:, 1]  
    Z_Sat = sol[:, 2]  

	# The `dataSet` variable is created as a 2D NumPy array containing the satellite's position components at each time step
    dataSet = np.array([X_Sat, Y_Sat, Z_Sat])  
    numDataPoints = len(t)  

	# The variables `N`, `phi`, and `theta` are defined to create a grid of points representing the Earth's surface in spherical coordinates
    N = 50  
    phi = np.linspace(0, 2 * np.pi, N)  
    theta = np.linspace(0, np.pi, N)  
    theta, phi = np.meshgrid(theta, phi)  
  
    r_Earth = 6378.14  

	# The `X_Earth`, `Y_Earth`, and `Z_Earth` arrays are computed using the spherical coordinates and the radius of the Earth
    X_Earth = r_Earth * np.cos(phi) * np.sin(theta)  
    Y_Earth = r_Earth * np.sin(phi) * np.sin(theta)  
    Z_Earth = r_Earth * np.cos(theta)  

	# The variables `x`, `y`, and `z` are defined to create a line animation effect
    x = np.sin(np.pi / 5 * t)  
    y = np.sin(np.pi / 3 * t)  
    z = np.linspace(0, 100, 100)  
  
    def animate_func(num):  

		# The 3D plot is cleared
        ax.clear()  

		# The satellite's trajectory up to the current frame is plotted
        ax.plot3D(dataSet[0, :num + 1], dataSet[1, :num + 1],  
                  dataSet[2, :num + 1], c='blue')  

		# The current position of the satellite is marked with a scatter plot
        ax.scatter(dataSet[0, num], dataSet[1, num], dataSet[2, num],  
                   c='blue', marker='o')  

		# The starting position of the satellite is plotted as a black marker
        ax.plot3D(dataSet[0, 0], dataSet[1, 0], dataSet[2, 0],  
                  c='black', marker='o')  
        ax.set_title('Trajectory \nTime = ' + str(np.round(t[num],  
                                                           decimals=2)) + ' s')  

		# The Earth's surface is plotted
        ax.plot_surface(X_Earth, Y_Earth, Z_Earth, color='blue', alpha=0.7)  
        ax.view_init(30, 145)  
        ax.set_xlabel('X [Km]')  
        ax.set_ylabel('Y [Km]')  
        ax.set_zlabel('Z [Km]')  
        xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),  
                           ax.get_zlim3d()]).T  
        XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])  
        ax.set_xlim3d(XYZlim)  
        ax.set_ylim3d(XYZlim)  
        ax.set_zlim3d(XYZlim * 3 / 4)  
  
    fig = plt.figure()  
    ax = plt.axes(projection='3d')  

	# The `animation.FuncAnimation()` function is used to animate the plot by repeatedly calling the `animate_func` function with a specified interval and number of frames
    line_ani = animation.FuncAnimation(fig, animate_func, interval=100,  
                                       frames=numDataPoints)  
    plt.show()
~~~
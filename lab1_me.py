import sympy as sp
from sympy import symbols, solve
from sympy.solvers import nonlinsolve
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.integrate import odeint

print("Task 1: Determination of the satellite's orbital elements")

def Task1 (x,y,z,Vx,Vy,Vz,M):
    mr = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    mv = sp.sqrt(Vx ** 2 + Vy ** 2 + Vz ** 2)
    r = sp.N(mr)
    v = sp.N(mv)
    sg = x * Vx + y * Vy + z * Vz
    gm = sp.acos(sg / (r * v))
    rp, vp = symbols('rp, vp', real=True)
    p = list(nonlinsolve([r * v * sp.sin(gm) - rp * vp, vp ** 2 - v ** 2 - 2 * M * (1 / rp - 1 / r)], [rp, vp]))
    A = p[0]
    B = p[1]
    if sp.N((r * v ** 2 / M - 1),1,chop=True)==0:
        p[0]=r
        p[1]=v
    elif A[0] > B[0] or A[0]<0:
        p = B
    else:
        p = A
    e = p[0] * p[1] ** 2 / M - 1
    e = sp.N(e, 5)

    nu = sp.atan((r * v ** 2 / M) * sp.sin(gm) * sp.cos(gm) / (r * v ** 2 / M * sp.sin(gm) ** 2 - 1))
    sg = v * r * sp.sin(gm)
    Cz = x * Vy - Vx * y
    Cy = z * Vx - x * Vz
    Cx = y * Vz - z * Vy
    ir = sp.acos(Cz / sg)
    i = ir * 180 / sp.pi
    i = sp.N(i, 8)
    Fx=Vy*Cz-Vz*Cy-M*x/r
    Fy=Vz*Cx-Vx*Cz-M*y/r
    Fz=Vx*Cy-Vy*Cx-M*z/r
    F=sp.sqrt(Fx**2+Fy**2+Fz**2)
    Enx=-Cy/(sp.sqrt(Cx**2+Cy**2))
    Eny=Cx/(sp.sqrt(Cx**2+Cy**2))
    omgr=sp.acos(((Enx*Fx)+(Eny*Fy))/F)
    omg = omgr * 180 / sp.pi
    omg = sp.N(omg,8)
    xr = sp.acos(-(Cy / sp.sqrt(Cx ** 2 + Cy ** 2)))
    x = xr * 180 / sp.pi
    x = sp.N(x,8)
    E0 = sp.acos((e + sp.cos(nu)) / (1 + e * sp.cos(nu)))
    EP = sp.acos((e + 1) / (1 + e))
    M0 = E0 - e * sp.sin(E0)
    MP = EP - e * sp.sin(EP)
    if sp.N(e,1)==1:
        Fp='Nan'
        Tp='Nan'
        a=0
        n=0
    else:
        a = p[0] / (1 - e)
        #print("a", a)
        #print("e", e)
        Fp = a * (1 - e ** 2)
        Fp = sp.N(Fp, 10)
        n = sp.sqrt(M / a ** 3)
        Tp = (MP - M0) / n
        if Tp<=0:
            Tp = sp.N(Tp,8)
        else:
            Tp = (M0 - MP) / n
            Tp = sp.N(Tp, 8)
    print('Focal Parameter =', Fp, "\nEccentricity =", e, "\nInclination =", i, '\nArgument of Periapsis =', omg,
          "\nLongitude of Ascending Node =", x)
    print('Time of passage of the periapsis t=', Tp, 'c', )
    return (Fp,e,ir,omgr,xr,M0,n,a)


def task2 (Fp,e,i,omg,x,N,M0,n,a,M):
    t1=N*60*60
    M1=M0+n*t1
    E=0
    Epr=-1
    while abs(E-Epr)>0.001:
        Epr=E
        E=e*sp.sin(E)+M1
    Nu=2*sp.atan(sp.sqrt((1+e)/(1-e))*sp.tan(E/2))
    r=(a*(1-e**2)/(1+e*sp.cos(Nu)))
    vt=(sp.sqrt(M/Fp)*(1+e*sp.cos(Nu)))
    vr=(sp.sqrt(M/Fp)*e*sp.sin(Nu))
    R=np.array([[r],[0],[0]])
    V=np.array([[vr],[vt],[0]])
    A=np.array([[(sp.cos(omg+Nu)*sp.cos(x)-sp.sin(omg+Nu)*sp.cos(i)*sp.sin(x)),(-sp.sin(omg+Nu)*sp.cos(x)-sp.cos(omg+Nu)*sp.cos(i)*sp.sin(x)),(sp.sin(i)*sp.sin(x))],\
                [(sp.cos(omg+Nu)*sp.sin(x))+sp.sin(omg+Nu)*sp.cos(i)*sp.cos(x),(-sp.sin(omg+Nu)*sp.sin(x)+sp.cos(omg+Nu)*sp.cos(i)*sp.cos(x)),(-sp.sin(i)*sp.cos(x))],\
                [(sp.sin(omg+Nu)*sp.sin(i)),(sp.cos(omg+Nu)*sp.sin(i)),(sp.cos(i))]])
    R=np.matmul(A,R)
    V=np.matmul(A,V)
    R=R.flatten().tolist()
    V=V.flatten().tolist()
    x,y,z=sp.N(R[0]),sp.N(R[1]),sp.N(R[2])
    Vx,Vy,Vz=sp.N(V[0]),sp.N(V[1]),sp.N(V[2])
    print('After ', N, ' hour(s)')
    print('Coordinates: х- ',x,' y- ',y,' z- ',z)
    print('Velocities: Vx- ',Vx,' Vy- ',Vy,' Vz- ',Vz)
    return E


def Topographic_coordinates(R, om_e, dt, lam0, R0):

    r = sp.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])
    dlam = sp.atan2(R[0], R[1]) - sp.atan2(R0[0], R0[1])
    lam = lam0 + dlam - om_e * dt

    while lam > sp.pi: lam -= 2 * sp.pi
    while lam < -sp.pi: lam += 2 * sp.pi

    fi = sp.asin(R[2] / r)
    return lam, fi


def task3 (Rr,lam0,a,M,e,i,omg,x,M0,n):
    print(a)
    T = 2 * sp.pi * sp.sqrt(a ** 3 / M)
    print("T", T)
    om_e=(360*sp.pi/180)/(23*60*60+56*60+4)
    t1=3*60
    dt=3*60
    r = sp.sqrt(Rr[0] * Rr[0] + Rr[1] * Rr[1] + Rr[2] * Rr[2])
    fi=sp.asin(Rr[2] / r)
    X,Y=[sp.N(lam0*180/sp.pi,4)],[sp.N(fi*180/sp.pi,4)]
    while t1<T*2:
        M1 = M0 + n * t1
        E = 0
        Epr = -1
        while abs(E - Epr) > 0.001:
            Epr = E
            E = e * sp.sin(E) + M1
        Nu = 2 * sp.atan(sp.sqrt((1 + e) / (1 - e)) * sp.tan(E / 2))
        r = (a * (1 - e ** 2) / (1 + e * sp.cos(Nu)))
        R = np.array([[r], [0], [0]])
        A = np.array([[(sp.cos(omg + Nu) * sp.cos(x) - sp.sin(omg + Nu) * sp.cos(i) * sp.sin(x)),
                       (-sp.sin(omg + Nu) * sp.cos(x) - sp.cos(omg + Nu) * sp.cos(i) * sp.sin(x)),
                       (sp.sin(i) * sp.sin(x))], \
                      [(sp.cos(omg + Nu) * sp.sin(x)) + sp.sin(omg + Nu) * sp.cos(i) * sp.cos(x),
                       (-sp.sin(omg + Nu) * sp.sin(x) + sp.cos(omg + Nu) * sp.cos(i) * sp.cos(x)),
                       (-sp.sin(i) * sp.cos(x))], \
                      [(sp.sin(omg + Nu) * sp.sin(i)), (sp.cos(omg + Nu) * sp.sin(i)), (sp.cos(i))]])
        R = np.matmul(A, R)
        R = R.flatten().tolist()
        R= [sp.N(R[0]), sp.N(R[1]), sp.N(R[2])]
        R0=Rr
        Xx,Yy=Topographic_coordinates(R, om_e, dt, lam0, R0)
        lam0 = Xx
        Xx=sp.N(Xx*180/sp.pi,4)
        Yy=sp.N(Yy*180/sp.pi,4)
        X.append(Xx)
        Y.append(Yy)
        Rr = R
        t1=t1+dt
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


def task3_3d (ZVx,ZVy,ZVz,ZX,ZY,ZZ,M,):
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

    X_0 = ZX
    Y_0 = ZY
    Z_0 = ZZ
    VX_0 = ZVx
    VY_0 = ZVy
    VZ_0 = ZVz
    state_0 = [X_0, Y_0, Z_0, VX_0, VY_0, VZ_0]
    t = np.linspace(0, 6 * 3600, 200)  # Simulates for a time period of 6
    # hours [s]
    sol = odeint(model_2BP, state_0, t)
    X_Sat = sol[:, 0]
    Y_Sat = sol[:, 1]
    Z_Sat = sol[:, 2]
    dataSet = np.array([X_Sat, Y_Sat, Z_Sat])
    numDataPoints = len(t)
    N = 50
    phi = np.linspace(0, 2 * np.pi, N)
    theta = np.linspace(0, np.pi, N)
    theta, phi = np.meshgrid(theta, phi)

    r_Earth = 6378.14
    X_Earth = r_Earth * np.cos(phi) * np.sin(theta)
    Y_Earth = r_Earth * np.sin(phi) * np.sin(theta)
    Z_Earth = r_Earth * np.cos(theta)
    x = np.sin(np.pi / 5 * t)
    y = np.sin(np.pi / 3 * t)
    z = np.linspace(0, 100, 100)

    def animate_func(num):
        ax.clear()
        ax.plot3D(dataSet[0, :num + 1], dataSet[1, :num + 1],
                  dataSet[2, :num + 1], c='blue')
        ax.scatter(dataSet[0, num], dataSet[1, num], dataSet[2, num],
                   c='blue', marker='o')
        ax.plot3D(dataSet[0, 0], dataSet[1, 0], dataSet[2, 0],
                  c='black', marker='o')
        ax.set_title('Trajectory \nTime = ' + str(np.round(t[num],
                                                           decimals=2)) + ' s')

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
    line_ani = animation.FuncAnimation(fig, animate_func, interval=100,
                                       frames=numDataPoints)
    plt.show()


x,y,z,Vx,Vy,Vz=-3200, 8200,5800,5,-2, 6
ZVx,ZVy,ZVz,ZX,ZY,ZZ= Vx,Vy,Vz,x,y,z
t=datetime(year=2024,month=12,day=3,hour=10,minute=0, second=0,  tzinfo=None)
N=1
M=398600
lam0=-4.8*sp.pi/180
Rr=[x,y,z]


Task1(x,y,z,Vx,Vy,Vz,M)
Fp,e,i,omg,x,M0,n,a=Task1(x,y,z,Vx,Vy,Vz,M)

print('')
print('Task 2: Using the orbit elements found in laboratory work No. 1, determine the coordinates and velocity components of the spacecraft in the geocentric equatorial coordinate system N hours after the moment to, where N is the number of your task.')
('')
task2(Fp,e,i,omg,x,N,M0,n,a,M)
print('\nTask 3: Construct a satellite track on two turns of its orbit')
print('')
X31,Y31=task3(Rr,lam0,a,M,e,i,omg,x,M0,n)
task3_3d(ZVx,ZVy,ZVz,ZX,ZY,ZZ,M,)


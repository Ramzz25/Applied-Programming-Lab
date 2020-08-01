#APL End semester
#Coded By : Ramanan S
#Date : 30/07/2020

#Importing required functions and libraries
from pylab import *
import math
import mpl_toolkits.mplot3d.axes3d as p3

#Calculating 'h', if required parameters are given
#Values are taken at SI units
#L->Inductance connected to circuit, A->Area of cross-section of plate
#w->Resonant frequency, Er->Dielectric constant
def calculate_h(Er,L,A,w,Ly):
    Eo = 8.85*(10**(-12))
    h = (((L*w*w*A*(Eo)) - Ly)*Er)/(1-Er)
    return h
print('Value of h for a given set of initial parameters :',calculate_h(2,0.015,1,10**6,0.2)*100,'cms')

#function to find the values of A and B using the error vector
def fitForError(errors, x):
    A = zeros((len(errors), 2))
    A[:, 0] = 1
    A[:, 1] = x
    return A, lstsq(A, log(errors), rcond=None)[0]

#Function to solve the Laplace equation taking the specified params
def solve_laplace(M,N,delta1,k,delta2,No,Er):
    Nx = M                      #size along x
    Ny = N                      #size along y
    Niter = No                  #No. of max iters
    phi = np.zeros((Ny,Nx))     #Defining potential matrix
    err = zeros(Niter)          #Defing array to record errors 
    itr_taken = Niter
    phi[0,:] = 1                #Giving +1v at the top
    for j in range(Niter):      #Calculating the potential values until convergence
        oldphi = phi.copy()
        
        phi[1:k, 1:-1] = 0.25*(phi[1:k, 0:-2]+phi[1:k, 2:]+phi[0:k-1, 1:-1]+phi[2:k+1, 1:-1])
        phi[k,1:-1] = (1/(1+Er))*( phi[k-1,1:-1] + Er*phi[k+1,1:-1] )
        phi[k+1:-1, 1:-1] = 0.25*(phi[k+1:-1, 0:-2]+phi[k+1:-1, 2:]+phi[k:-2, 1:-1]+phi[k+2:, 1:-1])
        
        err[j] = (abs(phi-oldphi)).max()
        
        
    M2, c2 = fitForError(err[500:], range(Niter)[500:]) # fit2 - taking iterations > 500
    A = (exp(c2[0]))
    B = c2[1]
    for j in range(500,Niter):
        er = (A*((exp(B*(j+0.5)))))/(-B)
        if(er<delta2):
            itr_taken=j+1
            break
            
    return phi,itr_taken,err

#Checking the working of the solve_laplace function by plotting it's output for a set of parameters
Ny = 400
Nx = 200
delta1 = 20/Ny
itr = 30000
k = 160
delta2 = (5.7)*(10**(-2))
Er = 2
phi,itr_taken,err = solve_laplace(Nx,Ny,delta1,k,delta2,itr,Er)
print("No. of iterations to reach the desired accuracy :",itr_taken)
#Plotting Error Vs No. of iterations
#We can observe error value converges
f = figure()
ax = f.add_subplot(111)
ax.plot(range(itr)[300::50], err[300::50], 'or')
title('Error Vs No. of iterations')
xlabel('No. of iters')
ylabel('Error')
show()

#Contour plot of potential of the tank for the above specified parameters
y = linspace(0, 20, Ny)  # range of y
x = linspace(0, 10, Nx)  # range of x
X, Y = meshgrid(x, -y)
fig = figure()
ax = fig.add_subplot(111)
plt1 = ax.contourf(X, Y, phi, cmap="RdBu_r")
title("Contour plot of Updated potential $\phi$")
xlabel("$x$")
ylabel("$y$")
ax = gca()
fig.colorbar(plt1, ax=ax, orientation='vertical')
show()

#Function to calculate chanrges on top plate and the side walls containg the fuid for a set of params
def cal_charge(Ny,Nx,delta1,itr,k,delta2,Er):
    
    #Finding the potential distribution for the given set of params
    phi,itr_taken,err = solve_laplace(Nx,Ny,delta1,k,delta2,itr,Er)
    
    #Distance between 2 cells
    delx = 10/Nx
    dely = 20/Ny
    
    #Summing up the electric field in the walls and bottom
    Q_fluid = 0
    E_nor1 = sum(phi[k+1:Ny,1])/delx
    E_nor2 = sum(phi[k+1:Ny,-2])/(delx)
    E_nor3 = sum(phi[Ny-2,1:Nx-1])/(dely)
    Q_fluid = E_nor1 + E_nor2 + E_nor3
    
    #summing up the electric field at the top
    Q_top = 0
    E_nor = sum(phi[0]-phi[1])/(dely)
    Q_top = E_nor
        
    return Q_fluid,Q_top

#Function to plot charge Vs h
def plot_charge():
    Qf,Qt = [], []
    
    #Calculating charge for different heights of liquid
    for i in range(1,10):
        ratio = i/10
        q1,q2 = cal_charge(200,100,20/Ny,10000,int(200-(ratio)*(200)),10**(-6),2)
        Qf.append(q1)
        Qt.append(q2)
        
    #Plotting Qtop, Qfluid Vs h (height of the fluid)
    f = figure(figsize=(8,10))
    ax = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    x = [i/10 for i in range(1,10)]
    ax.plot(x,Qf)
    ax2.plot(x,Qt)
    f.suptitle("Plots of Qtop vs h and Qfluid vs h")
    ax.set_xlabel('Ratio of h/Ly')
    ax.set_ylabel('Qfluid')
    ax.grid(True)
    ax2.set_xlabel('Ratio of h/Ly')
    ax2.set_ylabel('Qtop')
    ax2.grid(True)
    show()
    
plot_charge()

#Function to calculate field at the central point
def calc_field(Ny,Nx,delta1,itr,k,delta2,Er):
    
    #Calculating potential distribution
    phi,itr_taken,err = solve_laplace(Nx,Ny,delta1,k,delta2,itr,Er)
    
    delx = 10/Nx
    dely = 20/Ny
    
    #Calculating tangential and normal components of electric field in both medium
    Ey_air = (phi[k-1,Nx//2] - phi[k,Nx//2])/(dely)
    Ex_air = (phi[k,Nx//2] - phi[k,Nx//2 +1])/(delx)
    Ey_fluid = (phi[k,Nx//2] - phi[k+1,Nx//2])/(dely)
    Ex_fluid = (phi[k,Nx//2] - phi[k,Nx//2 +1])/(delx)
    
    print('Field at the centre')
    print('Normal Component of Electrin Field in \'air\' medium :',Ey_air,'V/cm')
    print('Normal Component of Electrin Field in \'Fluid\' medium :',Ey_fluid,'V/cm')
    print('Tangential Component of Electrin Field in \'air\' medium :',Ex_air,'V/cm')
    print('Tangential Component of Electrin Field in \'Fluid\' medium :',Ex_fluid,'V/cm')
    print()
    print('D in air :',Ey_air)
    print('D in Fluid :',Er*Ey_fluid)
    print()
    print('Ratio of fields (Ey_fluid/Ey_air) across the boundary')
    #Calculating field across the boundary
    Ey_a = (phi[k-1,1:-1] - phi[k,1:-1])/(dely)
    Ey_f = (phi[k,1:-1] - phi[k+1,1:-1])/(dely)
    print(Ey_f/Ey_a)
    
calc_field(200,100,20/Ny,18000,100,10**(-5),2)

#Function to verify snell's law
def verify_snell(Ny,Nx,delta1,itr,k,delta2,Er):
    
    phi,itr_taken,err = solve_laplace(Nx,Ny,delta1,k,delta2,itr,Er)
    
    delx = 10/Nx
    dely = 20/Ny
    
    #Calculating field at central point
    Ey_air = (phi[k-1,Nx//2] - phi[k,Nx//2])/(dely)
    Ex_air = (phi[k,Nx//2] - phi[k,Nx//2 +1])/(delx)
    Ey_fluid = (phi[k,Nx//2] - phi[k+1,Nx//2])/(dely)
    Ex_fluid = (phi[k,Nx//2] - phi[k,Nx//2 +1])/(delx)
    
    #Finding angle of incidence,refraction using tangential and normal components of field
    o1 = math.degrees(math.atan( 1/(Ey_air/Ex_air) ))
    o2 = math.degrees(math.atan( 1/(Ey_fluid/Ex_fluid) ))
    
    print('Angle of incidence (observed from air to fluid) :',o1,'deg')
    print('Angle of Refraction (observed from air to fluid) :',o2,'deg')
    
    #Calculating ratio of sine of angles
    u1 = math.sin(math.radians(o1))
    u2 = math.sin(math.radians(o2))
    
    #Ratio of refractive indices is propotional to sqrt(Er)
    print('sin(i)/sin(r) =',u1/u2)
    print('Ratio of Refractive indexes =',1/sqrt(Er))
    
verify_snell(200,100,20/Ny,18000,100,10**(-5),2)
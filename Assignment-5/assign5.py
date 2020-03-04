#APL Assignment - 5
#Done by - Ramanan S
#Date - 04/03/2020

#Importing required libraries
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import warnings
warnings.filterwarnings('ignore')

Nx=30 # size along x
Ny=30 # size along y
radius=0.35 # radius of central lead
Niter=2000 # number of iterations to perform

phi = np.zeros((Ny,Nx))        #initialising potential array
y = linspace(-0.5,0.5,Ny)      #defining x,y coordinates 
x = linspace(-0.5,0.5,Nx)
Y, X = meshgrid(y, x)

#initialising potential of the circle of radius 0.35 to 1V
ii = where(square(X) + square(Y) <= pow(radius, 2))
phi[ii] = 1.0

#plotting the potential contour of the copper plate at the initial state
fig, (ax1) = plt.subplots(nrows=1)
plt1 = ax1.contourf(Y, X, phi, cmap="RdBu_r")
title("Contour Plot of $\phi$")
xlabel("$x$")
ylabel("$y$")
ax = gca()
fig.colorbar(plt1, ax=ax, orientation='vertical')
show()

errors = zeros(Niter)  # initialise error array to zeros

for j in range(Niter):
    oldphi = phi.copy()  #making a copy of phi array to find the error
    
    #updating phi array ()
    phi[1:-1, 1:-1] = 0.25*(phi[1:-1, 0:-2]+phi[1:-1, 2:]+phi[0:-2, 1:-1]+phi[2:, 1:-1])

    phi[1:-1, 0] = phi[1:-1, 1]  # Left edge
    phi[1:-1, -1] = phi[1:-1, -2]  # right edge
    phi[0, :] = phi[1, :]  # Top edge

    # Assign 1 V to electrode region
    phi[ii] = 1.0
    
    # Appending errors for each iterations
    errors[j] = (abs(phi-oldphi)).max()

#plotting the error vs no. of iterations in both semilog and log-log axis
f1 = figure(figsize=(15,5))
ax = f1.add_subplot(121)
ax2 = f1.add_subplot(122)

ax.semilogy(range(Niter),errors,'r',label='err')
ax.title.set_text('semilog graph')
ax.grid()

ax2.loglog(range(Niter),errors,'r',label='err')
ax2.title.set_text('loglog graph')
ax2.grid()

show()

#plotting error after every 50th iteration
f = figure()
ax = f.add_subplot(111)
ax.plot(range(Niter)[300::50], errors[300::50], 'or')
show()

#plotting the 3-D surface plot of the potential at steady state
fig2=figure(figsize=(12,5)) # open a new figure
ax=p3.Axes3D(fig2) # Axes3D is the means to do a surface plot
title('The 3-D surface plot of the potential')
surf = ax.plot_surface(-Y, -X, phi, rstride=1, cstride=1, cmap=cm.jet)
fig2.colorbar(surf, orientation='vertical')
show()

#plotting the potential contour of the copper plate at steady state
fig3 = figure()
ax = fig3.add_subplot(111)
plt1 = ax.contourf(-Y, -X, phi, cmap="RdBu_r")
title("Contour plot of Updated potential $\phi$")
xlabel("$x$")
ylabel("$y$")
ax = gca()
fig.colorbar(plt1, ax=ax, orientation='vertical')
show()

#function to find the values of A and B using the error vector
def fitForError(errors, x):
    A = zeros((len(errors), 2))
    A[:, 0] = 1
    A[:, 1] = x
    return A, lstsq(A, log(errors), rcond=None)[0]

#function to compute error using log(A) and B values
def computeErrorFit(M, c):
    return exp(M.dot(c))

M1, c1 = fitForError(errors, range(Niter))  # fit1 - taking all iterations
M2, c2 = fitForError(errors[500:], range(Niter)[500:])  # fit2 - taking iterations > 500

print("Fit1 : A = %g , B = %g" % ((exp(c1[0]), c1[1])))
print("Fit2 : A = %g , B = %g" % ((exp(c2[0]), c2[1])))

print("The time Constant (1/B) all iterations: %g" % (abs(1/c1[1])))
print("The time Constant (1/B) for higher iterations (> 500) : %g" % (abs(1/c2[1])))

#calculating the error vector using the two fits
err_f1 = (exp(c1[0])*(exp(c1[1]*range(Niter))))
err_f2 = (exp(c2[0])*(exp(c2[1]*range(Niter))))

#plotting the calculated error using fits and original error
f4 = figure(figsize=(12,6))
ax = f4.add_subplot(111)

ax.semilogy(range(Niter),errors,'r',label='err')
ax.semilogy(range(Niter),err_f1,'c',label='err_f1')
ax.semilogy(range(Niter),err_f2,'--g',label='err_f2')

title('Errors')
ax.grid()
ax.legend()
show()

#initialising Jx, Jy arrays
Jx = zeros((Ny, Nx))
Jy = zeros((Ny, Nx))

#filling Jx, Jy arrays with conductivity = 1
Jx[1:-1, 1:-1] = 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])
Jy[1:-1, 1:-1] = 0.5*(phi[2:, 1:-1] - phi[0:-2, 1:-1])

#plotting J using quiver function
f5 = figure(figsize=(12,8))
ax = f5.add_subplot(111)

ax.scatter(x[ii[0]], y[ii[1]], color='r', label="V = 1 region")

ax.quiver(-Y, -X, Jx, Jy)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

ax.legend()
title("The Vector plot of the current flow")
show()
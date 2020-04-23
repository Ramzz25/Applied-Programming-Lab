#APL Assignment - 8
#Done BY : Ramanan S


# Importing required libraries.

from pylab import *

#defining a function transform, which plot the fft
def transform(func = sin, T = 8*pi, N = 512, lim = 5, ret = False, c = 1):
    t=linspace(-T/2,T/2,N+1);t=t[:-1]
    y=func(t)
    Y=fftshift(fft(ifftshift(y)))/N*c
    w=linspace(-N*(pi/T),N*(pi/T),N+1);w=w[:-1]
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=2)
    xlim([-lim,lim])
    ylabel(r"$|Y|$",size=16)
    grid(True)
    subplot(2,1,2)
    ii = where(abs(Y) > 1e-3)
    plot(w,angle(Y),'go',markersize= 3)
    plot(w[ii],angle(Y[ii]),'ro',markersize = 5)
    xlim([-lim,lim])
    ylabel(r"Phase of $Y$",size=16)
    xlabel(r"$\omega$",size=16)
    grid(True)
    show()
    if(ret):
        return Y

#plotting for sin^3(x) and cos^3(x)
transform(func = lambda x : (sin(x))**3, lim = 5)
transform(func = lambda x : (cos(x))**3, lim = 5)

#plotting for the function x(t) = cos(20t + 5cos(t))
transform(func = lambda x : cos(20*x + 5*cos(x)), lim = 40)

#defining and plotting for gaussian function
def gaussian_transform(x):
    return (1/sqrt(2*pi))*exp((-x**2)/2)

N = 128
T = 8 * pi
err = 1
tol = 1e-6
yprevious = 0
while(err > tol):
    print('For N = '+str(N)+' and T = '+str(T))
    y = transform(func = lambda x : exp(-x**2/2), T = T, N = N, ret = True, c = T/(2*pi))
    err = sum(abs(y[::2]-yprevious))
    yprevious = y
    N = 2*N
    T = 2*T

print('Preferred Values are :')
print('N = '+str(N))
print('T = '+str(T))
print('True Error = '+str(err))
print('The original transform is:')
w = linspace(-N*(pi/T),N*(pi/T),N+1);w=w[:-1]
y = gaussian_transform(w)
subplot(2,1,1)
plot(w,abs(y),lw=2)
xlim([-5,5])
ylabel(r"$|Y|$",size=16)
grid(True)
subplot(2,1,2)
grid()
ylabel("Phase of Y",size=16)
plot(w,angle(y),'go',markersize= 3)
show()